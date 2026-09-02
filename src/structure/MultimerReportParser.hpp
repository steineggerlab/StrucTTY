#pragma once

#include <string>
#include <vector>

struct MultimerHit {
    std::string qComplex;
    std::string tComplex;
    std::vector<std::string> qChains;
    std::vector<std::string> tChains;
    float qTMScore = -1.0f;
    float tTMScore = -1.0f;

    float U[9] = {1,0,0, 0,1,0, 0,0,1};
    float T[3] = {0,0,0};
    bool has_transform = false;

    float qComplexCov = -1.0f;
    float tComplexCov = -1.0f;
    std::string qChainTms;
    std::string tChainTms;
    std::string interfaceLddt;
    int assId = -1;
};

class MultimerReportParser {
public:
    bool load(const std::string& filepath);

    const std::vector<MultimerHit>& get_hits() const { return hits_; }
    int hit_count() const { return (int)hits_.size(); }

private:
    std::vector<MultimerHit> hits_;
};
