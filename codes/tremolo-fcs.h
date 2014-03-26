#ifndef TREMOLO_FCS_H
#define TREMOLO_FCS_H

#include <fcs.h>

#include "defs.h"
#include "parser.h"

struct FCSData;
struct Problem;

void sendFCSParametersToAll(struct FCSData *params);

void HandleFCSError(FCSResult fcs_result);

void StartFCS(struct Problem *P);

int ReadFCSCoulombRecord(struct Problem *P, FilePosType *filePos, parse_data *pd);

void InitFCS(struct Problem *P);

void UpdateFCSGrid(struct Problem *P);

void UpdateFCSStructFactors(struct Problem *P);

void FCSForce(struct Problem *P);

void CleanupFCS(struct Problem *P);

#endif /* !TREMOLO_FCS_H */
