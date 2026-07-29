#include <iostream>
#include <clocale>
#include <filesystem>
#include <memory>
#include <map>
#include <vector>
#include <string>
#include "structty.h"
#include "Terminal.hpp"
#include "Common.hpp"
#include "Screen.hpp"
#include "MSAParser.hpp"
#include "FoldseekParser.hpp"
#include "MultimerReportParser.hpp"
#include "FoldMasonParser.hpp"
#include "Benchmark.hpp"

namespace fs = std::filesystem;

namespace structty {

void run(const RunOptions& opts) {
    setlocale(LC_ALL, "");
    Terminal::enter_raw_mode();

    Terminal::Size term_sz = Terminal::get_size();
    int term_rows = term_sz.rows;
    int term_cols = term_sz.cols;

    Screen screen(term_cols, term_rows, opts.show_structure, opts.mode);

    Benchmark bm;
    const bool bench = opts.benchmark;
    using BenchClock = Benchmark::clock;
    BenchClock::time_point t_load0;

    if (bench) {
        fs::path in_file(opts.input_files[0]);
        bm.start(std::filesystem::current_path().string() + "/structty_bench_"
                 + current_timestamp() + "_" + in_file.stem().string() + ".csv");
        std::cout << "CWD = " << std::filesystem::current_path() << std::endl;
        screen.set_benchmark(&bm);
        t_load0 = Benchmark::clock::now();
    }

    screen.set_chainfile(opts.chains_file, (int)opts.input_files.size());

    // Query-from-DB (D3/D4) + multi-query nav (D7/D8/D9): when a query tmp DB is
    // provided, read the query structure(s) from the Foldseek query DB instead of
    // parsing the original CLI path. Handles folder/tar/gz query inputs that
    // gemmi's single-file reader cannot open. Hits are grouped by the .m8 query
    // column so ]/[ can navigate between queries.
    // Step 7 (D6): multimer `_report` 경로. 14컬럼 tsv 를 complex 단위로 파싱해
    // query/target complex 전체 체인을 로드하고 complex U/T 로 겹침. query/target DB 는
    // complex DB(체인별 엔트리). report_format 가 켜졌을 때만 진입(.m8 경로와 배타적).
    bool multimer_report = false;
    if (opts.report_format && !opts.foldseek_query_db.empty() &&
        !opts.foldseek_file.empty()) {
        MultimerReportParser mr_parser;
        if (mr_parser.load(opts.foldseek_file) && mr_parser.hit_count() > 0) {
            screen.set_multimer_report(mr_parser.get_hits(), opts.foldseek_query_db,
                                       opts.foldseek_db, opts.show_structure);
            multimer_report = true;
        } else {
            std::cerr << "Warning: failed to parse multimer report: "
                      << opts.foldseek_file << std::endl;
        }
    }

    bool query_from_db = false;
    std::vector<std::string> query_ids;                              // .m8 순서
    std::map<std::string, std::vector<FoldseekHit>> hits_by_query;   // query별 hit 그룹
    if (!multimer_report &&
        !opts.foldseek_query_db.empty() && !opts.foldseek_file.empty()) {
        FoldseekParser q_parser;
        if (q_parser.load(opts.foldseek_file) && q_parser.hit_count() > 0) {
            for (const FoldseekHit& h : q_parser.get_hits()) {
                auto it = hits_by_query.find(h.query);
                if (it == hits_by_query.end()) {
                    query_ids.push_back(h.query);
                    hits_by_query.emplace(h.query, std::vector<FoldseekHit>{ h });
                } else {
                    it->second.push_back(h);
                }
            }
            query_from_db = !query_ids.empty();
        }
    }

    if (multimer_report) {
        // Scene already set up by set_multimer_report() above (complex chains
        // loaded + superposed). Nothing more to do here.
    } else if (query_from_db) {
        // Screen owns the full per-query setup (load query from DB, normalize,
        // open target DB, load first hit). Query switching via ]/[ reuses the
        // same path. The workflow handoff passes only the query as input_files,
        // so plaintext targets are not loaded here.
        screen.set_query_nav(query_ids, hits_by_query, opts.foldseek_query_db,
                             opts.foldseek_db, opts.show_structure);
    } else {
    for (int i = 0; i < (int)opts.input_files.size(); i++) {
        screen.set_protein(opts.input_files[i], i, opts.show_structure);
    }
    screen.set_tmatrix();
    screen.normalize_proteins();
    screen.update_total_len_ca();

    // 기능 1: interface 모드일 때 inter-chain interface 계산 (threshold=8.0Å)
    if (opts.mode == "interface") {
        screen.compute_interface_all();
    }

    // 기능 4: aligned 모드 — foldseek/FoldMason 결과가 없는 경우 (nearest-neighbor fallback)
    if (opts.mode == "aligned" && opts.foldseek_file.empty() && opts.foldmason_file.empty()) {
        screen.compute_aligned_all();
    }

    // Foldseek DB 직접 읽기 모드 (--db 옵션)
    if (!opts.foldseek_db.empty()) {
        screen.open_foldseek_db(opts.foldseek_db);
    }

    // 기능 3: Foldseek hit 탐색 설정 (-fs 파일이 있을 때)
    if (!opts.foldseek_file.empty()) {
        FoldseekParser fs_nav_parser;
        if (fs_nav_parser.load(opts.foldseek_file) && fs_nav_parser.hit_count() > 0) {
            if ((int)opts.input_files.size() > 1) {
                // 작업 3-A: 다중 타겟 — m8 hit을 CLI target 파일명 기준으로 필터링
                std::vector<std::string> target_stems;
                for (int i = 1; i < (int)opts.input_files.size(); i++) {
                    fs::path p(opts.input_files[i]);
                    target_stems.push_back(p.stem().string());
                }

                const std::vector<FoldseekHit>& all_hits = fs_nav_parser.get_hits();
                std::vector<FoldseekHit> filtered_hits;
                for (const FoldseekHit& hit : all_hits) {
                    for (const std::string& stem : target_stems) {
                        if (hit.target.find(stem) != std::string::npos ||
                            stem.find(hit.target) != std::string::npos) {
                            filtered_hits.push_back(hit);
                            break;
                        }
                    }
                }
                screen.set_foldseek_hits(filtered_hits);
                screen.set_fs_db_path(opts.foldseek_db_path);
                screen.prepare_foldseek_db(filtered_hits);

                // 작업 3-B: 각 target protein에 매칭 hit의 transform 적용
                for (int ti = 1; ti < (int)opts.input_files.size(); ti++) {
                    fs::path tp(opts.input_files[ti]);
                    std::string tstem = tp.stem().string();

                    for (const FoldseekHit& hit : filtered_hits) {
                        if (hit.target.find(tstem) != std::string::npos ||
                            tstem.find(hit.target) != std::string::npos) {
                            screen.apply_hit_transform(ti, hit);
                            break;
                        }
                    }
                }
            } else {
                // 단일 입력: 기존 동작 유지
                screen.set_foldseek_hits(fs_nav_parser.get_hits());
                screen.set_fs_db_path(opts.foldseek_db_path);
                screen.prepare_foldseek_db(fs_nav_parser.get_hits());
                screen.load_next_hit(+1);
            }
        }
    }

    // 기능 5: conservation 모드일 때 MSA 파일 로드 및 conservation score 계산
    if (opts.mode == "conservation" && !opts.msa_file.empty()) {
        MSAParser msa_parser;
        if (msa_parser.load(opts.msa_file)) {
            std::vector<float> scores = msa_parser.compute_conservation();
            screen.apply_msa_conservation(0, scores);
        } else {
            std::cerr << "Warning: Failed to load MSA file: " << opts.msa_file << std::endl;
        }
    }

    // 기능 8: FoldMason MSA 기반 superposition + conservation
    if (!opts.foldmason_file.empty()) {
        auto fm_parser = std::make_unique<FoldMasonParser>();
        const std::string& fm_path = opts.foldmason_file;
        bool fm_loaded = false;

        std::string ext;
        {
            size_t dot = fm_path.rfind('.');
            if (dot != std::string::npos) ext = fm_path.substr(dot + 1);
            for (auto& c : ext) c = (char)std::tolower((unsigned char)c);
        }

        if (ext == "json") {
            fm_loaded = fm_parser->load_json(fm_path);
        } else {
            fm_loaded = fm_parser->load_fasta(fm_path);
        }

        if (!fm_loaded || fm_parser->entry_count() == 0) {
            std::cerr << "Warning: Failed to load FoldMason file: " << fm_path << std::endl;
        } else {
            // entry 매칭: 파일명(확장자 제외) 기준으로 탐색, 실패 시 0/1 순서 가정
            int fm_query_idx = 0;
            int fm_target_idx = 1;
            const auto& entries = fm_parser->get_entries();
            if ((int)opts.input_files.size() >= 1) {
                fs::path qp(opts.input_files[0]);
                std::string qstem = qp.stem().string();
                for (int ei = 0; ei < fm_parser->entry_count(); ei++) {
                    if (entries[ei].name.find(qstem) != std::string::npos) {
                        fm_query_idx = ei; break;
                    }
                }
            }
            if ((int)opts.input_files.size() >= 2) {
                fs::path tp(opts.input_files[1]);
                std::string tstem = tp.stem().string();
                for (int ei = 0; ei < fm_parser->entry_count(); ei++) {
                    if (ei == fm_query_idx) continue;
                    if (entries[ei].name.find(tstem) != std::string::npos) {
                        fm_target_idx = ei; break;
                    }
                }
            }

            int total_entries = fm_parser->entry_count();

            // conservation 색상 (단일/다중 구조 모두)
            if (opts.mode == "conservation") {
                std::vector<float> entropy = fm_parser->compute_column_entropy(false);
                auto col_map = fm_parser->build_query_col_map(fm_query_idx);
                std::vector<float> mapped_scores(col_map.size(), 0.0f);
                for (int ri = 0; ri < (int)col_map.size(); ri++) {
                    int col = col_map[ri];
                    if (col < (int)entropy.size()) mapped_scores[ri] = entropy[col];
                }
                screen.apply_msa_conservation(0, mapped_scores);
            }

            // 패널 정보 구성
            FoldMasonInfo fm_info;
            fm_info.valid = true;
            fm_info.entry_count = total_entries;

            // 두 구조 superposition (두 번째 PDB 있을 때만)
            bool do_superposition = ((int)opts.input_files.size() >= 2 && total_entries >= 2);
            if (do_superposition) {
                fm_info.align_method = (opts.mode == "aligned") ? "msa-col" : "-";
            } else {
                fm_info.align_method = "-";
            }

            screen.set_foldmason(std::move(fm_parser));
            screen.set_foldmason_panel_info(fm_info);

            if (do_superposition) {
                screen.apply_foldmason_superposition(0, 1, fm_query_idx, fm_target_idx);
            }
        }
    }
    } // end else (non-query-from-DB scene setup)

    if (bench) {
        auto t_load1 = Benchmark::clock::now();
        bm.log("load", -1, Benchmark::ms_since(t_load0, t_load1));

        const std::vector<int> script = {
            'X','Y','Z','A','D','W','S','R','F'
        };

        const int warmup = 200;
        const int events = 2000;

        // Warmup run (not measured)
        bool old_enabled = bm.enabled;
        bm.enabled = false;
        for (int i = 0; i < warmup; i++) {
            screen.draw_screen(opts.no_panel);
            screen.handle_input(script[i % script.size()]);
        }
        bm.enabled = old_enabled;

        // Measured run
        bm.log("bench_begin", -1, 0.0);
        for (int i = 0; i < events; i++) {
            screen.draw_screen(opts.no_panel);
            screen.handle_input(script[i % script.size()]);
        }
        bm.log("bench_end", -1, 0.0);
    } else {
        bool running = true;
        bool needs_redraw = true;
        while (running) {
            if (needs_redraw) {
                screen.draw_screen(opts.no_panel);
            }
            // KEY_MOUSE 이벤트 시 needs_redraw=false: 패널 부분 갱신만 수행
            // 키보드 이벤트 시 needs_redraw=true: 전체 재렌더링
            running = screen.handle_input(needs_redraw);
        }
    }

    Terminal::exit_raw_mode();
}

} // namespace structty
