#include "InputProbe.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace input_probe {

namespace {

std::string to_lower(const std::string& s) {
    std::string out = s;
    std::transform(out.begin(), out.end(), out.begin(),
                   [](unsigned char c) { return (char)std::tolower(c); });
    return out;
}

bool ends_with(const std::string& s, const std::string& suffix) {
    if (s.size() < suffix.size()) return false;
    return s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool has_structure_extension(const std::string& path) {
    std::string p = to_lower(path);
    if (ends_with(p, ".gz")) {
        p.erase(p.size() - 3);
    }
    return ends_with(p, ".pdb") || ends_with(p, ".cif") || ends_with(p, ".ent");
}

bool exists_file(const std::string& path) {
    std::error_code ec;
    return fs::is_regular_file(path, ec);
}

bool is_directory(const std::string& path) {
    std::error_code ec;
    return fs::is_directory(path, ec);
}

std::string first_data_line(const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs.is_open()) return std::string();

    std::string line;
    while (std::getline(ifs, line)) {
        const size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) continue;
        if (line[first] == '#') continue;
        return line;
    }
    return std::string();
}

InputKind kind_from_columns(int ncols) {
    if (ncols == 14) return InputKind::ResultReport;
    if (ncols == 12 || ncols == 17 || ncols >= 21) return InputKind::ResultM8;
    return InputKind::Unknown;
}

std::string representative_file(const std::string& dir) {
    std::error_code ec;
    std::vector<std::string> names;
    for (const auto& entry : fs::directory_iterator(dir, ec)) {
        if (ec) return std::string();
        std::error_code file_ec;
        if (!entry.is_regular_file(file_ec)) continue;
        names.push_back(entry.path().string());
    }
    if (names.empty()) return std::string();
    return *std::min_element(names.begin(), names.end());
}

}

int tsv_column_count(const std::string& path) {
    const std::string line = first_data_line(path);
    if (line.empty()) return 0;

    int fields = 1;
    for (const char c : line) {
        if (c == '\t') ++fields;
    }
    return fields;
}

bool db_has_ca(const std::string& db_path) {
    return exists_file(db_path + "_ca.dbtype");
}

InputKind probe(const std::string& path) {
    if (path.empty()) return InputKind::Unknown;

    if (exists_file(path + ".dbtype")) return InputKind::FoldseekDB;

    if (is_directory(path)) {
        const std::string sample = representative_file(path);
        if (sample.empty()) return InputKind::Unknown;
        if (has_structure_extension(sample)) return InputKind::StructureDir;
        const std::string line = first_data_line(sample);
        if (!line.empty() && line[0] == '>') return InputKind::SequenceFasta;
        return InputKind::Unknown;
    }

    if (!exists_file(path)) return InputKind::Unknown;

    if (has_structure_extension(path)) return InputKind::Structure;

    const std::string line = first_data_line(path);
    if (line.empty()) return InputKind::Unknown;

    if (line[0] == '>') return InputKind::SequenceFasta;

    return kind_from_columns(tsv_column_count(path));
}

const char* kind_name(const InputKind kind) {
    switch (kind) {
        case InputKind::Structure:     return "structure file";
        case InputKind::StructureDir:  return "structure directory";
        case InputKind::FoldseekDB:    return "Foldseek DB";
        case InputKind::ResultM8:      return "Foldseek m8";
        case InputKind::ResultReport:  return "Foldseek multimer report";
        case InputKind::SequenceFasta: return "sequence FASTA";
        case InputKind::Unknown:       break;
    }
    return "unknown";
}

namespace {

bool check_has_coordinates(const std::string& path, const std::string& role) {
    const InputKind kind = probe(path);

    if (kind == InputKind::SequenceFasta) {
        std::cerr << "Error: " << role << " '" << path << "' is a sequence FASTA.\n"
                  << "       StrucTTY needs 3D coordinates: pass PDB/mmCIF files, a directory\n"
                  << "       of them, or a Foldseek DB built from structures.\n"
                  << "       (foldseek createdb --prostt5-model predicts 3Di but writes no _ca)"
                  << std::endl;
        return false;
    }
    if (kind == InputKind::FoldseekDB && !db_has_ca(path)) {
        std::cerr << "Error: Foldseek DB '" << path << "' has no C-alpha coordinates ("
                  << path << "_ca is missing).\n"
                  << "       It was built from sequences, or with --index-exclude 2."
                  << std::endl;
        return false;
    }
    return true;
}

}

bool validate_inputs(const std::vector<std::string>& query,
                     const std::string& target, const std::string& result) {
    if (query.empty()) {
        std::cerr << "Error: Need input file dir" << std::endl;
        return false;
    }

    if (!target.empty() && result.empty()) {
        std::cerr << "Error: -fst / --foldseek-target given without -fsr / --foldseek-result.\n"
                  << "       Both must be given together." << std::endl;
        return false;
    }
    if (target.empty() && !result.empty()) {
        std::cerr << "Error: -fsr / --foldseek-result given without -fst / --foldseek-target.\n"
                  << "       Both must be given together. Use -fst auto to download hit structures."
                  << std::endl;
        return false;
    }

    for (const std::string& q : query) {
        if (!check_has_coordinates(q, "query")) {
            return false;
        }
    }
    if (!target.empty() && target != "auto" && !check_has_coordinates(target, "-fst target")) {
        return false;
    }

    if (result.empty()) {
        return true;
    }

    const InputKind result_kind = probe(result);
    if (result_kind != InputKind::ResultM8 && result_kind != InputKind::ResultReport) {
        const int ncols = tsv_column_count(result);
        std::cerr << "Error: -fsr '" << result << "' is not a Foldseek result file.\n"
                  << "       Expected tab-separated 12/17/21/29 columns (m8) or 14 columns"
                  << " (multimer _report),\n       ";
        if (ncols == 0) {
            std::cerr << "but the file could not be read or has no data rows.";
        } else {
            std::cerr << "but found " << ncols << " columns.";
        }
        std::cerr << std::endl;
        return false;
    }

    return true;
}

}
