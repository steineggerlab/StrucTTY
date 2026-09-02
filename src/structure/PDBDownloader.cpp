#include "PDBDownloader.hpp"

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>
#include <sys/stat.h>
#include <cerrno>

static bool file_exists_nonempty(const std::string& path) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0) return false;
    return st.st_size > 0;
}

static bool mkdir_p(const std::string& path) {
    struct stat st;
    if (stat(path.c_str(), &st) == 0) return true;
#if defined(_WIN32)
    return _mkdir(path.c_str()) == 0 || errno == EEXIST;
#else
    return mkdir(path.c_str(), 0755) == 0 || errno == EEXIST;
#endif
}

static std::string get_home_dir() {
    const char* home = std::getenv("HOME");
    if (home) return std::string(home);
    return "/tmp";
}

static std::string cache_root() {
    return get_home_dir() + "/.cache/structty/pdb";
}

static void ensure_cache_dir() {
    std::string home = get_home_dir();
    mkdir_p(home + "/.cache");
    mkdir_p(home + "/.cache/structty");
    mkdir_p(home + "/.cache/structty/pdb");
}

static std::string to_lower(const std::string& s) {
    std::string r = s;
    for (char& c : r) c = (char)std::tolower((unsigned char)c);
    return r;
}

static bool starts_with(const std::string& s, const std::string& prefix) {
    return s.size() >= prefix.size() && s.substr(0, prefix.size()) == prefix;
}

static bool contains(const std::string& s, const std::string& sub) {
    return s.find(sub) != std::string::npos;
}

static bool is_alnum_upper(char c) {
    return (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9');
}

static bool is_alnum_lower_or_digit(char c) {
    return (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9');
}

std::string PDBDownloader::extract_target_id(const std::string& raw_target) {
    size_t pos = raw_target.find(' ');
    if (pos == std::string::npos) return raw_target;
    return raw_target.substr(0, pos);
}

DBType PDBDownloader::detect_db_type(const std::string& target_id) {
    if (target_id.empty()) return DBType::Unknown;

    if (starts_with(target_id, "AF-") && contains(target_id, "_TED")) {
        size_t ted_pos = target_id.rfind("_TED");
        if (ted_pos != std::string::npos) {
            std::string after = target_id.substr(ted_pos + 4);
            bool all_digit = !after.empty();
            for (char c : after) if (c < '0' || c > '9') { all_digit = false; break; }
            if (all_digit) return DBType::TED;
        }
    }

    if (contains(target_id, "_unrelaxed_rank_001_alphafold2")) {
        return DBType::BFVD_Local;
    }

    if (starts_with(target_id, "GMGC10.")) {
        return DBType::GMGCL;
    }

    if (starts_with(target_id, "BFD") && target_id.size() > 3) {
        size_t i = 3;
        while (i < target_id.size() && target_id[i] >= '0' && target_id[i] <= '9') i++;
        if (i > 3 && i < target_id.size() && target_id[i] == '_') {
            return DBType::BFMD;
        }
    }

    if (starts_with(target_id, "AF-")) {
        size_t dash1 = target_id.find('-', 3);
        if (dash1 != std::string::npos) {
            size_t dash2 = target_id.find('-', dash1 + 1);
            if (dash2 != std::string::npos) {
                std::string frag = target_id.substr(dash1 + 1, dash2 - dash1 - 1);
                if (!frag.empty() && frag[0] == 'F' &&
                    target_id.find("model_v", dash2) != std::string::npos) {
                    return DBType::AlphaFold;
                }
            }
        }
    }

    if (starts_with(target_id, "MGYP") && target_id.size() == 16) {
        bool all_digit = true;
        for (size_t i = 4; i < 16; i++) {
            if (target_id[i] < '0' || target_id[i] > '9') { all_digit = false; break; }
        }
        if (all_digit) return DBType::ESMAtlas30;
    }

    if (target_id.size() >= 6 &&
        target_id[0] >= '0' && target_id[0] <= '9' &&
        is_alnum_lower_or_digit(target_id[1]) &&
        is_alnum_lower_or_digit(target_id[2]) &&
        is_alnum_lower_or_digit(target_id[3]) &&
        target_id[4] == '-' &&
        contains(target_id, ".cif.gz_")) {
        return DBType::PDB;
    }

    if (target_id.size() >= 6 &&
        target_id[0] >= '0' && target_id[0] <= '9' &&
        is_alnum_lower_or_digit(target_id[1]) &&
        is_alnum_lower_or_digit(target_id[2]) &&
        is_alnum_lower_or_digit(target_id[3]) &&
        target_id[4] == '_' &&
        target_id.size() >= 6) {
        return DBType::PDB;
    }

    if (target_id.size() == 6 &&
        target_id[0] >= '0' && target_id[0] <= '9' &&
        is_alnum_lower_or_digit(target_id[1]) &&
        is_alnum_lower_or_digit(target_id[2]) &&
        is_alnum_lower_or_digit(target_id[3]) &&
        std::isalpha((unsigned char)target_id[4]) &&
        target_id[5] >= '0' && target_id[5] <= '9') {
        return DBType::CATH50;
    }

    {
        std::string base = target_id;
        size_t underscore = target_id.rfind('_');
        if (underscore != std::string::npos) {
            std::string suffix = target_id.substr(underscore + 1);
            bool all_digit = !suffix.empty();
            for (char c : suffix) if (c < '0' || c > '9') { all_digit = false; break; }
            if (all_digit) base = target_id.substr(0, underscore);
        }
        if (base.size() >= 6 && base.size() <= 10) {
            bool all_upper = true;
            for (char c : base) if (!is_alnum_upper(c)) { all_upper = false; break; }
            if (all_upper) return DBType::BFVD_Official;
        }
    }

    return DBType::Unknown;
}

std::string PDBDownloader::extract_chain(const std::string& target_id, DBType db_type) {
    if (db_type != DBType::PDB) return "-";

    size_t cifgz_pos = target_id.find(".cif.gz_");
    if (cifgz_pos != std::string::npos) {
        return target_id.substr(cifgz_pos + 8);
    }

    size_t pos = target_id.find('_');
    if (pos == std::string::npos) return "-";

    std::string after = target_id.substr(pos + 1);
    size_t last_us = after.rfind('_');
    if (last_us != std::string::npos) {
        after = after.substr(last_us + 1);
    }
    return after.empty() ? "-" : after;
}

std::string PDBDownloader::extract_pdb_id(const std::string& target_id) {
    size_t dash_pos = target_id.find('-');
    if (dash_pos == 4) return target_id.substr(0, 4);

    size_t pos = target_id.find('_');
    if (pos == std::string::npos) return target_id.substr(0, 4);
    return target_id.substr(0, std::min(pos, (size_t)4));
}

std::string PDBDownloader::extract_uniprot_id(const std::string& target_id, DBType db_type) {
    if (db_type == DBType::BFVD_Local) {
        size_t pos = target_id.find('_');
        if (pos != std::string::npos) return target_id.substr(0, pos);
        return target_id;
    }
    if (db_type == DBType::BFVD_Official) {
        size_t underscore = target_id.rfind('_');
        if (underscore != std::string::npos) {
            std::string suffix = target_id.substr(underscore + 1);
            bool all_digit = !suffix.empty();
            for (char c : suffix) if (c < '0' || c > '9') { all_digit = false; break; }
            if (all_digit) return target_id.substr(0, underscore);
        }
        return target_id;
    }
    return target_id;
}

std::string PDBDownloader::get_download_url(const std::string& target_id, DBType db_type) {
    switch (db_type) {
        case DBType::PDB: {
            std::string pdb_id = extract_pdb_id(target_id);
            return "https://files.rcsb.org/download/" + pdb_id + ".cif";
        }
        case DBType::AlphaFold: {
            return "https://alphafold.ebi.ac.uk/files/" + target_id + ".cif";
        }
        case DBType::ESMAtlas30: {
            return "https://api.esmatlas.com/fetchPredictedStructure/" + target_id;
        }
        case DBType::CATH50: {
            return "https://www.cathdb.info/version/v4_3_0/api/rest/id/" + target_id + ".pdb";
        }
        case DBType::BFVD_Local:
        case DBType::BFVD_Official: {
            std::string uniprot = extract_uniprot_id(target_id, db_type);
            return "https://bfvd.steineggerlab.workers.dev/pdb/" + uniprot + ".pdb";
        }
        case DBType::TED: {
            return "https://ted.cathdb.info/api/v1/files/" + target_id + ".pdb";
        }
        case DBType::GMGCL:
        case DBType::BFMD:
        case DBType::Unknown:
        default:
            return "";
    }
}

std::string PDBDownloader::get_cache_path(const std::string& target_id, DBType db_type) {
    ensure_cache_dir();
    std::string root = cache_root();

    switch (db_type) {
        case DBType::PDB: {
            std::string pdb_id = extract_pdb_id(target_id);
            return root + "/" + pdb_id + ".cif";
        }
        case DBType::AlphaFold:
            return root + "/" + target_id + ".cif";
        case DBType::ESMAtlas30:
            return root + "/" + target_id + ".pdb";
        case DBType::CATH50:
            return root + "/" + target_id + ".pdb";
        case DBType::BFVD_Local:
        case DBType::BFVD_Official: {
            std::string uniprot = extract_uniprot_id(target_id, db_type);
            return root + "/" + uniprot + ".pdb";
        }
        case DBType::TED:
            return root + "/" + target_id + ".pdb";
        default:
            return root + "/" + target_id + ".pdb";
    }
}

static std::string find_by_name(const std::string& name, const std::string& db_path) {
    const std::string base       = db_path + "/" + name;
    const std::string base_lower = db_path + "/" + to_lower(name);

    const std::vector<std::string> candidates = {
        base + ".pdb",       base + ".cif",       base + ".ent",
        base + ".pdb.gz",    base + ".cif.gz",    base + ".ent.gz",
        base,
        base_lower + ".pdb", base_lower + ".cif", base_lower + ".ent",
        base_lower + ".pdb.gz", base_lower + ".cif.gz", base_lower + ".ent.gz",
        base_lower,
    };

    for (const std::string& cand : candidates) {
        if (file_exists_nonempty(cand)) return cand;
    }
    return "";
}

std::string PDBDownloader::find_in_db_path(const std::string& target_id,
                                            const std::string& db_path) {
    if (db_path.empty()) return "";

    const std::string direct = find_by_name(target_id, db_path);
    if (!direct.empty()) return direct;

    size_t cut = target_id.rfind('_');
    while (cut != std::string::npos && cut > 0) {
        const std::string found = find_by_name(target_id.substr(0, cut), db_path);
        if (!found.empty()) return found;
        if (cut == 0) break;
        cut = target_id.rfind('_', cut - 1);
    }
    return "";
}

bool PDBDownloader::download_file(const std::string& url, const std::string& dest_path,
                                  const std::string& extra_header) {
    if (url.empty() || dest_path.empty()) return false;

    std::string cmd;
    if (extra_header.empty()) {
        cmd = "curl -f -s -L -o \"" + dest_path + "\" \"" + url + "\" 2>/dev/null"
              " || wget -q -O \"" + dest_path + "\" \"" + url + "\" 2>/dev/null";
    } else {
        cmd = "curl -f -s -L -H \"" + extra_header + "\" -o \"" + dest_path + "\" \"" + url + "\" 2>/dev/null"
              " || wget -q --header=\"" + extra_header + "\" -O \"" + dest_path + "\" \"" + url + "\" 2>/dev/null";
    }

    FILE* p = popen(cmd.c_str(), "r");
    if (!p) return false;
    pclose(p);

    return file_exists_nonempty(dest_path);
}

std::string PDBDownloader::get_no_url_message(DBType db_type,
                                               const std::string& target_id) {
    switch (db_type) {
        case DBType::GMGCL:
            return "GMGCL: no download URL available. Web-server only DB.";
        case DBType::BFMD:
            return "BFMD: no download URL. Use --db-path.";
        default:
            return "File not found: " + target_id;
    }
}

std::string PDBDownloader::resolve_target_file(const std::string& target_id,
                                               const std::string& db_path,
                                               std::string& status_msg) {
    status_msg.clear();
    DBType db_type = detect_db_type(target_id);

    if (!db_path.empty()) {
        std::string found = find_in_db_path(target_id, db_path);
        if (!found.empty()) return found;
    }

    std::string cache_path = get_cache_path(target_id, db_type);
    if (file_exists_nonempty(cache_path)) return cache_path;

    std::string url = get_download_url(target_id, db_type);
    if (!url.empty()) {
        std::string header;
        if (db_type == DBType::TED) {
            header = "accept: application/octet-stream";
        }
        bool ok = download_file(url, cache_path, header);
        if (ok) return cache_path;
        if (db_type == DBType::ESMAtlas30) {
            status_msg = "ESMAtlas API rate limit, retry later";
        } else {
            status_msg = "Download failed: " + target_id;
        }
        return "";
    }

    status_msg = get_no_url_message(db_type, target_id);
    return "";
}
