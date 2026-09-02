#pragma once

#include <string>

enum class DBType {
    PDB,
    AlphaFold,
    ESMAtlas30,
    CATH50,
    BFVD_Local,
    BFVD_Official,
    GMGCL,
    TED,
    BFMD,
    Unknown
};

class PDBDownloader {
public:
    static std::string extract_target_id(const std::string& raw_target);

    static DBType detect_db_type(const std::string& target_id);

    static std::string extract_chain(const std::string& target_id, DBType db_type);

    static std::string extract_pdb_id(const std::string& target_id);

    static std::string extract_uniprot_id(const std::string& target_id, DBType db_type);

    static std::string get_download_url(const std::string& target_id, DBType db_type);

    static std::string get_cache_path(const std::string& target_id, DBType db_type);

    static std::string find_in_db_path(const std::string& target_id,
                                       const std::string& db_path);

    static bool download_file(const std::string& url, const std::string& dest_path,
                              const std::string& extra_header = "");

    static std::string resolve_target_file(const std::string& target_id,
                                           const std::string& db_path,
                                           std::string& status_msg);

    static std::string get_no_url_message(DBType db_type, const std::string& target_id);
};
