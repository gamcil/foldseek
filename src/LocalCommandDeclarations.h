#ifndef LOCALCOMMANDDECLARATIONS_H
#define LOCALCOMMANDDECLARATIONS_H

#include "Command.h"

extern int structuresearch(int argc, const char** argv, const Command &command);
extern int structcreatedb(int argc, const char** argv, const Command &command);
extern int structureindex(int argc, const char** argv, const Command &command);
extern int structurecluster(int argc, const char** argv, const Command &command);
extern int easystructuresearch(int argc, const char** argv, const Command &command);
extern int easystructurecluster(int argc, const char** argv, const Command &command);
extern int easymsa(int argc, const char** argv, const Command &command);
extern int tmalign(int argc, const char** argv, const Command &command);
extern int aln2tmscore(int argc, const char** argv, const Command &command);
extern int structurealign(int argc, const char** argv, const Command &command);
extern int samplemulambda(int argc, const char** argv, const Command &command);
extern int structureconvertalis(int argc, const char** argv, const Command &command);
extern int structureto3didescriptor(int argc, const char** argv, const Command &command);
extern int structurerbh(int argc, const char** argv, const Command &command);
extern int structureeasyrbh(int argc, const char** argv, const Command &command);
extern int structureungappedalign(int argc, const char** argv, const Command &command);
extern int structuremsa(int argc, const char** argv, const Command &command);
extern int msa2lddt(int argc, const char** argv, const Command &command);
#endif
