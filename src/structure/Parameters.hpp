#pragma once
#include <iostream>
#include <cstring>
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <vector>

using namespace std;
namespace fs = std::filesystem;

class Parameters{
    private:
        bool show_structure = false;
        bool no_panel = false;
        bool predict_structure = false;
        bool arg_okay = true;
        bool benchmark_mode = false;
        vector<string> in_file;
        vector<string> chains;
        string utmatrix = "";
        string chainfile = "";
        string mode = "protein";
        string depthcharacter = "@#%*^-.";
        string msa_file = "";
        string foldseek_file = "";
        string db_path = "";
    public:
        Parameters(int argc, char* argv[]);

        void print_args();

        bool is_valid_number(const std::string& str, int min, int max);

        // get, set
        vector<string>& get_in_file(){
            return in_file;
        }
        string get_in_file(int idx){
            if (idx < in_file.size()){
                return in_file[idx];
            }
            return "";
        }
        string get_chainfile(){
            return chainfile;
        }
        string get_utmatrix(){
            return utmatrix;
        }
        string get_mode(){
            return mode;
        }
        string get_depthcharacter(){
            return depthcharacter;
        }
        bool get_show_structure(){
            return show_structure;
        }
        bool get_no_panel(){
            return no_panel;
        }
        bool get_benchmark_mode(){
            return benchmark_mode;
        }
        bool check_arg_okay(){
            return arg_okay;
        }
        string get_msa_file(){
            return msa_file;
        }
        string get_foldseek_file(){
            return foldseek_file;
        }
        string get_db_path(){
            return db_path;
        }
};