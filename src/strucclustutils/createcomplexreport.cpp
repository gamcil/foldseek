#include "Util.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "IndexReader.h"
#include "FileUtil.h"
#include "TranslateNucl.h"
#include "MemoryMapped.h"
#include "Coordinate16.h"
#include "createcomplexreport.h"
#define ZSTD_STATIC_LINKING_ONLY
#include <zstd.h>
#include "LDDT.h"
#include "CalcProbTP.h"
#include <map>
#include <map>

#ifdef OPENMP
#include <omp.h>
#endif

typedef std::pair<std::string, std::string> compNameChainName_t;
typedef std::pair<unsigned int, std::string> ComplexAlignmentKey_t;

void getComplexNameChainName(std::string &chainName, compNameChainName_t &compAndChainName) {
    size_t pos = chainName.rfind('_');
    std::string comp = chainName.substr(0, pos);
    std::string chain = chainName.substr(pos + 1);
    compAndChainName = {comp, chain};
}

void getResult (std::vector<std::string> &qChainVector, std::vector<std::string> &tChainVector, std::vector<ComplexResult> &complexResVec, float qTMScore, float tTMScore, int assId) {
    char buffer[1024];
    std::string result;
    std::string qComplexName;
    std::string tComplexName;
    std::string qChainString;
    std::string tChainString;
    compNameChainName_t compAndChainName;
    getComplexNameChainName(qChainVector[0], compAndChainName);
    qComplexName = compAndChainName.first;
    qChainString = compAndChainName.second;
    getComplexNameChainName(tChainVector[0], compAndChainName);
    tComplexName = compAndChainName.first;
    tChainString = compAndChainName.second;
    for (size_t qChainId = 1; qChainId < qChainVector.size(); qChainId++) {
        getComplexNameChainName(qChainVector[qChainId], compAndChainName);
        qChainString += ',' + compAndChainName.second;
    }
    for (size_t tChainId = 1; tChainId < tChainVector.size(); tChainId++) {
        getComplexNameChainName(tChainVector[tChainId], compAndChainName);
        tChainString += ',' + compAndChainName.second;
    }
    int count = snprintf(buffer,sizeof(buffer),"%s\t%s\t%s\t%s\t%1.5f\t%1.5f\t%d\n", qComplexName.c_str(), tComplexName.c_str(), qChainString.c_str(), tChainString.c_str(), qTMScore, tTMScore, assId);
    result.append(buffer, count);
    complexResVec.emplace_back(ComplexResult(assId, result));
}

struct ComplexAlignment {
    ComplexAlignment(){};
    ComplexAlignment(std::string qChain, std::string tChain, double qTMscore, double tTMscore) : qTMScore(qTMscore),tTMScore(tTMscore){
        qChainVector = {qChain};
        tChainVector = {tChain};
    };
    std::vector<std::string> qChainVector;
    std::vector<std::string> tChainVector;
    double qTMScore;
    double tTMScore;
};


int createcomplexreport(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    const bool sameDB = par.db1.compare(par.db2) == 0 ? true : false;
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    int dbaccessMode = (DBReader<unsigned int>::USE_INDEX);
    std::map<unsigned int, unsigned int> qKeyToSet;
    std::map<unsigned int, unsigned int> tKeyToSet;
    IndexReader qDbr(par.db1, par.threads,  IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
    IndexReader qDbrHeader(par.db1, par.threads, IndexReader::SRC_HEADERS , (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
    IndexReader *tDbrHeader;
    if(sameDB){
        tDbrHeader = &qDbrHeader;
    } else{
        tDbrHeader = new IndexReader(par.db2, par.threads, IndexReader::SRC_HEADERS, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
    }

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    size_t localThreads = 1;
#ifdef OPENMP
    localThreads = std::max(std::min((size_t)par.threads, alnDbr.getSize()), (size_t)1);
#endif

    const bool shouldCompress = par.dbOut == true && par.compressed == true;
    const int dbType = par.dbOut == true ? Parameters::DBTYPE_GENERIC_DB : Parameters::DBTYPE_OMIT_FILE;
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), localThreads, shouldCompress, dbType);
    resultWriter.open();
    const bool isDb = par.dbOut;
    TranslateNucl translateNucl(static_cast<TranslateNucl::GenCode>(par.translationTable));
    Debug::Progress progress(alnDbr.getSize());
    std::vector<ComplexResult> complexResVec;
    Matcher::result_t res;
    auto complexDataHandler = ComplexDataHandler();
    std::map<ComplexAlignmentKey_t, ComplexAlignment> complexAlignmentsWithAssId;

#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#pragma omp  for schedule(dynamic, 10)
        for (size_t i = 0; i < alnDbr.getSize(); i++) {
            progress.updateProgress();
            const unsigned int queryKey = alnDbr.getDbKey(i);
            size_t qHeaderId = qDbrHeader.sequenceReader->getId(queryKey);
            const char *qHeader = qDbrHeader.sequenceReader->getData(qHeaderId, thread_idx);
            compNameChainName_t qCompAndChainName;
            std::string queryId = Util::parseFastaHeader(qHeader);
            getComplexNameChainName(queryId, qCompAndChainName);
            char *data = alnDbr.getData(i, thread_idx);
            while (*data != '\0') {
                bool isValid = parseScoreComplexResult(data, res, complexDataHandler);
                if (!isValid) {
                    std::cout << "error message";
                }
                data = Util::skipLine(data);
                size_t tHeaderId = tDbrHeader->sequenceReader->getId(res.dbKey);
                const char *tHeader = tDbrHeader->sequenceReader->getData(tHeaderId, thread_idx);
                std::string targetId = Util::parseFastaHeader(tHeader);
                unsigned int assId = complexDataHandler.assId;
                auto key = ComplexAlignmentKey_t(assId, qCompAndChainName.first);
                if (complexAlignmentsWithAssId.find(key) == complexAlignmentsWithAssId.end()){
                    complexAlignmentsWithAssId.insert({key, ComplexAlignment(queryId, targetId, complexDataHandler.qTmScore, complexDataHandler.tTmScore)});
                } else {
                    complexAlignmentsWithAssId[key].qChainVector.emplace_back(queryId);
                    complexAlignmentsWithAssId[key].tChainVector.emplace_back(targetId);
                }
            } // while end
        } // for end
    }
    std::map<ComplexAlignmentKey_t, ComplexAlignment>::iterator iter;
    for (iter = complexAlignmentsWithAssId.begin(); iter != complexAlignmentsWithAssId.end(); iter++) {
        getResult(iter->second.qChainVector, iter->second.tChainVector, complexResVec, iter->second.qTMScore, iter->second.tTMScore, iter->first.first);
    }
    SORT_SERIAL(complexResVec.begin(), complexResVec.end(), compareComplexResult);
    for (size_t i=0; i < complexResVec.size(); i++) {
        resultWriter.writeData(complexResVec[i].result.c_str(), complexResVec[i].result.length(), 0, localThreads - 1, false, false);
    }
    resultWriter.close(true);
    if (isDb == false) {
        FileUtil::remove(par.db4Index.c_str());
    }
    alnDbr.close();
    if (sameDB == false) {
        delete tDbrHeader;
    }
    return EXIT_SUCCESS;
}