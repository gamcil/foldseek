//
// Created by Woosub Kim on 2023/06/20.
//

#ifndef FOLDSEEK_CREATECOMPLEXREPORT_H
#define FOLDSEEK_CREATECOMPLEXREPORT_H
#include "Matcher.h"

const unsigned int NOT_AVAILABLE_CHAIN_KEY = 4294967295;
const double MAX_ASSIGNED_CHAIN_RATIO = 1.0;
const double TOO_SMALL_MEAN = 1.0;
const double TOO_SMALL_CV = 0.1;
const double FILTERED_OUT = 0.0;
const unsigned int UNCLUSTERED = 0;
const unsigned int CLUSTERED = 1;
const unsigned int MAX_RECURSIVE_NUM = 1000;
const unsigned int MIN_PTS = 2;

typedef std::pair<std::string, std::string> compNameChainName_t;
typedef std::map<unsigned int, unsigned int> chainKeyToComplexId_t;
typedef std::map<unsigned int, std::vector<unsigned int>> complexIdToChainKeys_t;
typedef std::vector<std::vector<unsigned int>> cluster_t;
typedef std::map<std::pair<unsigned int, unsigned int>, double> distMap_t;
typedef std::string resultToWrite_t;
typedef std::pair<unsigned int, resultToWrite_t> resultToWriteWithKey_t;

struct ScoreComplexResult {
    ScoreComplexResult() {}
    ScoreComplexResult(unsigned int assId, resultToWrite_t resultToWrite) : assId(assId), resultToWrite(resultToWrite) {}
    unsigned int assId;
    resultToWrite_t resultToWrite;
};

static bool compareComplexResult(const ScoreComplexResult &first, const ScoreComplexResult &second) {
    if (first.assId < second.assId) {
        return true;
    }
    if (first.assId > second.assId) {
        return false;
    }
    return false;
}

struct ComplexDataHandler {
    ComplexDataHandler(): assId(UINT_MAX), qTmScore(0.0f), tTmScore(0.0f) {}
    ComplexDataHandler(bool isValid): assId(UINT_MAX), qTmScore(0.0f), tTmScore(0.0f), isValid(isValid) {}
    ComplexDataHandler(unsigned int assId, double qTmScore, double tTmScore, std::string uString, std::string tString, bool isValid)
            : assId(assId), qTmScore(qTmScore), tTmScore(tTmScore), uString(uString), tString(tString), isValid(isValid) {}
    unsigned int assId;
    double qTmScore;
    double tTmScore;
    std::string uString;
    std::string tString;
    bool isValid;
};


static void getKeyToIdMapIdToKeysMapIdVec(
        const std::string &file,
        std::map<unsigned int, unsigned int> &chainKeyToComplexIdLookup,
        std::map<unsigned int, std::vector<unsigned int>> &complexIdToChainKeysLookup,
        std::vector<unsigned int> &complexIdVec
) {
    if (file.length() == 0) return;
    MemoryMapped lookupDB(file, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    char *data = (char *) lookupDB.getData();
    const char *entry[255];
    int prevComplexId =  -1;
    while (*data != '\0') {
        const size_t columns = Util::getWordsOfLine(data, entry, 255);
        if (columns < 3) {
            Debug(Debug::WARNING) << "Not enough columns in lookup file " << file << "\n";
            continue;
        }
        auto chainKey = Util::fast_atoi<int>(entry[0]);
        auto complexId = Util::fast_atoi<int>(entry[2]);
        chainKeyToComplexIdLookup.emplace(chainKey, complexId);
        if (complexId != prevComplexId) {
            complexIdToChainKeysLookup.emplace(complexId, std::vector<unsigned int>());
            complexIdVec.emplace_back(complexId);
            prevComplexId = complexId;
        }
        complexIdToChainKeysLookup.at(complexId).emplace_back(chainKey);
        data = Util::skipLine(data);
    }
    lookupDB.close();
}

static ComplexDataHandler parseScoreComplexResult(const char *data, Matcher::result_t &res) {
    const char *entry[255];
    size_t columns = Util::getWordsOfLine(data, entry, 255);
    if (columns!=16) {
        return ComplexDataHandler(false);
    }
    char key[255];
    ptrdiff_t keySize =  (entry[1] - data);
    strncpy(key, data, keySize);
    key[keySize] = '\0';
    unsigned int dbKey = Util::fast_atoi<unsigned int>(key);
    int score = Util::fast_atoi<int>(entry[1]);
    double seqId = strtod(entry[2],NULL);
    double eval = strtod(entry[3],NULL);
    int qStartPos =  Util::fast_atoi<int>(entry[4]);
    int qEndPos = Util::fast_atoi<int>(entry[5]);
    int qLen = Util::fast_atoi<int>(entry[6]);
    int dbStartPos = Util::fast_atoi<int>(entry[7]);
    int dbEndPos = Util::fast_atoi<int>(entry[8]);
    int dbLen = Util::fast_atoi<int>(entry[9]);
    int adjustQstart = (qStartPos==-1)? 0 : qStartPos;
    int adjustDBstart = (dbStartPos==-1)? 0 : dbStartPos;
    auto backtrace = std::string(entry[10], entry[11] - entry[10]);
    double qCov = SmithWaterman::computeCov(adjustQstart, qEndPos, qLen);
    double dbCov = SmithWaterman::computeCov(adjustDBstart, dbEndPos, dbLen);
    size_t alnLength = Matcher::computeAlnLength(adjustQstart, qEndPos, adjustDBstart, dbEndPos);
    double qTmScore = std::stod(entry[11]);
    double tTmScore = std::stod(entry[12]);
    std::string uString = std::string(entry[13], entry[14] - entry[13]-1);
    std::string tString = std::string(entry[14], entry[15] - entry[14]-1);
    unsigned int assId = Util::fast_atoi<unsigned int>(entry[15]);
    res = Matcher::result_t(dbKey, score, qCov, dbCov, seqId, eval, alnLength, qStartPos, qEndPos, qLen, dbStartPos, dbEndPos, dbLen, -1, -1, -1, -1, backtrace);
    return ComplexDataHandler(assId, qTmScore, tTmScore, uString, tString, true);
}

#endif //FOLDSEEK_CREATECOMPLEXREPORT_H
