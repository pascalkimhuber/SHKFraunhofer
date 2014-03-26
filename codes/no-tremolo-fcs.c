#include "no-tremolo-fcs.h"
#include "parser.h"

int ReadFCSCoulombRecord(struct Problem *P, FilePosType *filePos, parse_data *pd)
{
  (void)P;
  (void)pd;
  parse_error(ERRPARSE, filePos, "ReadFCSCoulombRecord: FCS not found/not implemented sequentially");
}

void StartFCS(struct Problem *P)
{
  (void)P;
}

void InitFCS(struct Problem *P)
{
  (void)P;
}

void UpdateFCSGrid(struct Problem *P)
{
  (void)P;
}

void UpdateFCSStructFactors(struct Problem *P)
{
  (void)P;
}

void FCSForce(struct Problem *P)
{
  (void)P;
}

void CleanupFCS(struct Problem *P)
{
  (void)P;
}
