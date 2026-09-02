#pragma once

#include <string>
#include <vector>

namespace input_probe {

enum class InputKind {
    Structure,
    StructureDir,
    FoldseekDB,
    ResultM8,
    ResultReport,
    SequenceFasta,
    Unknown
};

InputKind probe(const std::string& path);

bool db_has_ca(const std::string& db_path);

int tsv_column_count(const std::string& path);

const char* kind_name(InputKind kind);

bool validate_inputs(const std::vector<std::string>& query,
                     const std::string& target, const std::string& result);

}
