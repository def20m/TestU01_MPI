#include <unistd.h>
#include <chrono>
#include <cstring>
#include <string>
extern "C"
{
#include "TestU01.h"
#include "util.h"
#include "unif01.h"
#include "bbattery.h"
#include "swrite.h"
}

using std::string;

#define LEN 120
#define NAMELEN 30
#define NDIM 200

#define THOUSAND 1000
#define MILLION (THOUSAND * THOUSAND)
#define BILLION (THOUSAND * MILLION)

#define SMALLCRUSH_NUM 10
#define CRUSH_NUM 96
#define BIGCRUSH_NUM 106

std::string selectSmallCrushTest(unif01_Gen *gen, int testNum);
std::string selectCrushTest(unif01_Gen *gen, int testNum);
std::string selectBigCrushTest(unif01_Gen *gen, int testNum);

std::string GetPVal_Walk_nu(long N, swalk_Res *res);
std::string smarsa_BirthdaySpacings_nu(unif01_Gen *gen, sres_Poisson *res, long N, long n, int r, long d, int t, int p);
std::string sknuth_Collision_nu(unif01_Gen *gen, sknuth_Res2 *res, long N, long n, int r, long d, int t);
std::string sknuth_Gap_nu(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, double Alpha, double Beta);
std::string sknuth_SimpPoker_nu(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, int d, int k);
std::string sknuth_CouponCollector_nu(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, int d);
std::string sknuth_MaxOft_nu_Mean_Mean(unif01_Gen *gen, sknuth_Res1 *res, long N, long n, int r, int d, int t);
std::string sknuth_MaxOft_nu_Sum_AD(unif01_Gen *gen, sknuth_Res1 *res, long N, long n, int r, int d, int t);
std::string svaria_WeightDistrib_nu(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, long k, double alpha, double beta);
std::string swalk_RandomWalk1_nu(unif01_Gen *gen, swalk_Res *res, long N, long n, int r, int s, long L0, long L1, long pN);
std::string GetPVal_CPairs_nu(long N, snpair_Res *res);
std::string smarsa_SerialOver_nu(unif01_Gen *gen, sres_Basic *res, long N, long n, int r, long d, int t);
std::string smarsa_CollisionOver_nu(unif01_Gen *gen, smarsa_Res *res, long N, long n, int r, long d, int t);
std::string snpair_ClosePairs_GetPVal_CPairs_nu(unif01_Gen *gen, snpair_Res *res, long N, long n, int r, int t, int p, int m, long pN);
std::string snpair_ClosePairsBitMatch_nu(unif01_Gen *gen, snpair_Res *res, long N, long n, int r, int t);
std::string sknuth_Run_nu_Mean(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, lebool Up);
std::string sknuth_Run_nu_Sum(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, lebool Up);
std::string sknuth_Permutation_nu(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, int t);
std::string sknuth_CollisionPermut_nu(unif01_Gen *gen, sknuth_Res2 *res, long N, long n, int r, int t);
std::string svaria_SampleProd_nu_Mean(unif01_Gen *gen, sres_Basic *res, long N, long n, int r, int t);
std::string svaria_SampleProd_nu_AD(unif01_Gen *gen, sres_Basic *res, long N, long n, int r, int t);
std::string svaria_SampleMean_nu(unif01_Gen *gen, sres_Basic *res, long N, long n, int r);
std::string svaria_SampleCorr_nu(unif01_Gen *gen, sres_Basic *res, long N, long n, int r, int k);
std::string svaria_AppearanceSpacings_nu(unif01_Gen *gen, sres_Basic *res, long N, long Q, long K, int r, int s, int L);
std::string svaria_SumCollector_nu(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, double g);
std::string smarsa_MatrixRank_nu_Mean(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, int s, int L, int k);
std::string smarsa_MatrixRank_nu_Sum(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, int s, int L, int k);
std::string smarsa_Savir2_nu_Mean(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, long m, int t);
std::string smarsa_Savir2_nu_Sum(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, long m, int t);
std::string smarsa_GCD_nu_Mean(unif01_Gen *gen, smarsa_Res2 *res, long N, long n, int r, int s);
std::string smarsa_GCD_nu_Sum(unif01_Gen *gen, smarsa_Res2 *res, long N, long n, int r, int s);
std::string scomp_LinearComp_nu(unif01_Gen *gen, scomp_Res *res, long N, long n, int r, int s);
std::string scomp_LempelZiv_nu(unif01_Gen *gen, sres_Basic *res, long N, int k, int r, int s);
std::string sspectral_Fourier3_nu(unif01_Gen *gen, sspectral_Res *res, long N, int k, int r, int s);
std::string sstring_LongestHeadRun_nu(unif01_Gen *gen, sstring_Res2 *res, long N, long n, int r, int s, long L);
std::string sstring_PeriodsInStrings_nu_Mean(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, int s);
std::string sstring_PeriodsInStrings_nu_Sum(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, int s);
std::string sstring_HammingWeight2_nu(unif01_Gen *gen, sres_Basic *res, long N, long n, int r, int s, long L);
std::string sstring_HammingCorr_nu(unif01_Gen *gen, sstring_Res *res, long N, long n, int r, int s, int L);
std::string sstring_HammingIndep_nu_Mean(unif01_Gen *gen, sstring_Res *res, long N, long n, int r, int s, int L, int d);
std::string sstring_HammingIndep_nu_Sum(unif01_Gen *gen, sstring_Res *res, long N, long n, int r, int s, int L, int d);
std::string sstring_Run_nu(unif01_Gen *gen, sstring_Res3 *res, long N, long n, int r, int s);
std::string sstring_AutoCor_nu(unif01_Gen *gen, sres_Basic *res, long N, long n, int r, int s, int d);

std::string smarsa_BirthdaySpacings_nu(unif01_Gen *gen, sres_Poisson *res, long N, long n, int r, long d, int t, int p)
{
    smarsa_BirthdaySpacings(gen, res, N, n, r, d, t, p);
    std::string pVal = std::to_string(res->pVal2);
    return pVal;
}

std::string sknuth_Collision_nu(unif01_Gen *gen, sknuth_Res2 *res, long N, long n, int r, long d, int t)
{
    sknuth_Collision(gen, res, N, n, r, d, t);
    std::string pVal = std::to_string(res->Pois->pVal2);
    return pVal;
}

std::string sknuth_Gap_nu(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, double Alpha, double Beta)
{
    sknuth_Gap(gen, res, N, n, r, Alpha, Beta);
    std::string pVal = std::to_string(res->pVal2[gofw_Mean]);
    return pVal;
}

std::string sknuth_SimpPoker_nu(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, int d, int k)
{
    sknuth_SimpPoker(gen, res, N, n, r, d, k);
    std::string pVal = std::to_string(res->pVal2[gofw_Mean]);
    return pVal;
}

std::string sknuth_CouponCollector_nu(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, int d)
{
    sknuth_CouponCollector(gen, res, N, n, r, d);
    std::string pVal = std::to_string(res->pVal2[gofw_Mean]);
    return pVal;
}

std::string sknuth_MaxOft_nu_Mean_Mean(unif01_Gen *gen, sknuth_Res1 *res, long N, long n, int r, int d, int t)
{
    sknuth_MaxOft(gen, res, N, n, r, d, t);
    std::string pVal = std::to_string(res->Chi->pVal2[gofw_Mean]) + ",";
    pVal += std::to_string(res->Bas->pVal2[gofw_Mean]);
    return pVal;
}

std::string sknuth_MaxOft_nu_Sum_AD(unif01_Gen *gen, sknuth_Res1 *res, long N, long n, int r, int d, int t)
{

    sknuth_MaxOft(gen, res, N, n, r, d, t);
    std::string pVal = std::to_string(res->Chi->pVal2[gofw_Sum]) + ",";
    pVal += std::to_string(res->Bas->pVal2[gofw_AD]);
    return pVal;
}

std::string svaria_WeightDistrib_nu(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, long k, double alpha, double beta)
{
    svaria_WeightDistrib(gen, res, N, n, r, k, alpha, beta);
    std::string pVal = std::to_string(res->pVal2[gofw_Mean]);
    return pVal;
}

std::string smarsa_MatrixRank_nu_Mean(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, int s, int L, int k)
{
    smarsa_MatrixRank(gen, res, N, n, r, s, L, k);
    std::string pVal = std::to_string(res->pVal2[gofw_Mean]);
    return pVal;
}

std::string smarsa_MatrixRank_nu_Sum(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, int s, int L, int k)
{
    smarsa_MatrixRank(gen, res, N, n, r, s, L, k);
    std::string pVal = std::to_string(res->pVal2[gofw_Sum]);
    return pVal;
}

std::string sstring_HammingIndep_nu_Mean(unif01_Gen *gen, sstring_Res *res, long N, long n, int r, int s, int L, int d)
{
    sstring_HammingIndep(gen, res, N, n, r, s, L, d);
    std::string pVal = std::to_string(res->Bas->pVal2[gofw_Mean]);
    return pVal;
}

std::string sstring_HammingIndep_nu_Sum(unif01_Gen *gen, sstring_Res *res, long N, long n, int r, int s, int L, int d)
{
    sstring_HammingIndep(gen, res, N, n, r, s, L, d);
    std::string pVal = std::to_string(res->Bas->pVal2[gofw_Sum]);
    return pVal;
}

std::string swalk_RandomWalk1_nu(unif01_Gen *gen, swalk_Res *res, long N, long n, int r, int s, long L0, long L1, long pN)
{
    swalk_RandomWalk1(gen, res, N, n, r, s, L0, L1);
    std::string pVal = GetPVal_Walk_nu(pN, res);
    return pVal;
}

std::string smarsa_SerialOver_nu(unif01_Gen *gen, sres_Basic *res, long N, long n, int r, long d, int t)
{
    smarsa_SerialOver(gen, res, N, n, r, d, t);
    std::string pVal = std::to_string(res->pVal2[gofw_Mean]);
    return pVal;
}

std::string smarsa_CollisionOver_nu(unif01_Gen *gen, smarsa_Res *res, long N, long n, int r, long d, int t)
{
    smarsa_CollisionOver(gen, res, N, n, r, d, t);
    std::string pVal = std::to_string(res->Pois->pVal2);
    return pVal;
}

std::string snpair_ClosePairs_GetPVal_CPairs_nu(unif01_Gen *gen, snpair_Res *res, long N, long n, int r, int t, int p, int m, long pN)
{
    snpair_ClosePairs(gen, res, N, n, r, t, p, m);
    std::string pVal = GetPVal_CPairs_nu(pN, res);
    return pVal;
}

std::string snpair_ClosePairsBitMatch_nu(unif01_Gen *gen, snpair_Res *res, long N, long n, int r, int t)
{

    snpair_ClosePairsBitMatch(gen, res, N, n, r, t);
    std::string pVal = std::to_string(res->pVal[snpair_BM]);
    return pVal;
}

std::string sknuth_Run_nu_Sum(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, lebool Up)
{

    sknuth_Run(gen, res, N, n, r, Up);
    std::string pVal = std::to_string(res->pVal2[gofw_Sum]);
    return pVal;
}

std::string sknuth_Run_nu_Mean(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, lebool Up)
{

    sknuth_Run(gen, res, N, n, r, Up);
    std::string pVal = std::to_string(res->pVal2[gofw_Mean]);
    return pVal;
}

std::string sknuth_Permutation_nu(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, int t)
{

    sknuth_Permutation(gen, res, N, n, r, t);
    std::string pVal = std::to_string(res->pVal2[gofw_Mean]);
    return pVal;
}

std::string sknuth_CollisionPermut_nu(unif01_Gen *gen, sknuth_Res2 *res, long N, long n, int r, int t)
{

    sknuth_CollisionPermut(gen, res, N, n, r, t);
    std::string pVal = std::to_string(res->Pois->pVal2);
    return pVal;
}

std::string svaria_SampleProd_nu_Mean(unif01_Gen *gen, sres_Basic *res, long N, long n, int r, int t)
{

    svaria_SampleProd(gen, res, N, n, r, t);
    std::string pVal = std::to_string(res->pVal2[gofw_Mean]);
    return pVal;
}

std::string svaria_SampleProd_nu_AD(unif01_Gen *gen, sres_Basic *res, long N, long n, int r, int t)
{

    svaria_SampleProd(gen, res, N, n, r, t);
    std::string pVal = std::to_string(res->pVal2[gofw_AD]);
    return pVal;
}

std::string svaria_SampleMean_nu(unif01_Gen *gen, sres_Basic *res, long N, long n, int r)
{

    svaria_SampleMean(gen, res, N, n, r);
    std::string pVal = std::to_string(res->pVal2[gofw_AD]);
    return pVal;
}

std::string svaria_SampleCorr_nu(unif01_Gen *gen, sres_Basic *res, long N, long n, int r, int k)
{

    svaria_SampleCorr(gen, res, N, n, r, k);
    std::string pVal = std::to_string(res->pVal2[gofw_Mean]);
    return pVal;
}

std::string svaria_AppearanceSpacings_nu(unif01_Gen *gen, sres_Basic *res, long N, long Q, long K, int r, int s, int L)
{

    svaria_AppearanceSpacings(gen, res, N, Q, K, r, s, L);
    std::string pVal = std::to_string(res->pVal2[gofw_Mean]);
    return pVal;
}

std::string svaria_SumCollector_nu(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, double g)
{
    svaria_SumCollector(gen, res, N, n, r, g);
    std::string pVal = std::to_string(res->pVal2[gofw_Mean]);
    return pVal;
}

std::string smarsa_Savir2_nu_Mean(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, long m, int t)
{
    smarsa_Savir2(gen, res, N, n, r, m, t);
    std::string pVal = std::to_string(res->pVal2[gofw_Mean]);
    return pVal;
}

std::string smarsa_Savir2_nu_Sum(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, long m, int t)
{
    smarsa_Savir2(gen, res, N, n, r, m, t);
    std::string pVal = std::to_string(res->pVal2[gofw_Sum]);
    return pVal;
}

std::string smarsa_GCD_nu_Mean(unif01_Gen *gen, smarsa_Res2 *res, long N, long n, int r, int s)
{
    smarsa_GCD(gen, res, N, n, r, s);
    std::string pVal = std::to_string(res->GCD->pVal2[gofw_Mean]);
    return pVal;
}

std::string smarsa_GCD_nu_Sum(unif01_Gen *gen, smarsa_Res2 *res, long N, long n, int r, int s)
{
    smarsa_GCD(gen, res, N, n, r, s);
    std::string pVal = std::to_string(res->GCD->pVal2[gofw_Sum]);
    return pVal;
}

std::string scomp_LinearComp_nu(unif01_Gen *gen, scomp_Res *res, long N, long n, int r, int s)
{
    scomp_LinearComp(gen, res, N, n, r, s);
    std::string pVal = std::to_string(res->JumpNum->pVal2[gofw_Mean]) + ",";
    pVal += std::to_string(res->JumpSize->pVal2[gofw_Mean]);
    return pVal;
}

std::string scomp_LempelZiv_nu(unif01_Gen *gen, sres_Basic *res, long N, int k, int r, int s)
{

    scomp_LempelZiv(gen, res, N, k, r, s);
    std::string pVal = std::to_string(res->pVal2[gofw_Sum]);
    return pVal;
}

std::string sspectral_Fourier3_nu(unif01_Gen *gen, sspectral_Res *res, long N, int k, int r, int s)
{

    sspectral_Fourier3(gen, res, N, k, r, s);
    std::string pVal = std::to_string(res->Bas->pVal2[gofw_AD]);
    return pVal;
}

std::string sstring_LongestHeadRun_nu(unif01_Gen *gen, sstring_Res2 *res, long N, long n, int r, int s, long L)
{

    sstring_LongestHeadRun(gen, res, N, n, r, s, L);
    std::string pVal = std::to_string(res->Chi->pVal2[gofw_Mean]) + ",";
    pVal += std::to_string(res->Disc->pVal2);
    return pVal;
}

std::string sstring_PeriodsInStrings_nu_Mean(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, int s)
{

    sstring_PeriodsInStrings(gen, res, N, n, r, s);
    std::string pVal = std::to_string(res->pVal2[gofw_Mean]);
    return pVal;
}

std::string sstring_PeriodsInStrings_nu_Sum(unif01_Gen *gen, sres_Chi2 *res, long N, long n, int r, int s)
{

    sstring_PeriodsInStrings(gen, res, N, n, r, s);
    std::string pVal = std::to_string(res->pVal2[gofw_Sum]);
    return pVal;
}

std::string sstring_HammingWeight2_nu(unif01_Gen *gen, sres_Basic *res, long N, long n, int r, int s, long L)
{

    sstring_HammingWeight2(gen, res, N, n, r, s, L);
    std::string pVal = std::to_string(res->pVal2[gofw_Sum]);
    return pVal;
}

std::string sstring_HammingCorr_nu(unif01_Gen *gen, sstring_Res *res, long N, long n, int r, int s, int L)
{

    sstring_HammingCorr(gen, res, N, n, r, s, L);
    std::string pVal = std::to_string(res->Bas->pVal2[gofw_Mean]);
    return pVal;
}

std::string sstring_Run_nu(unif01_Gen *gen, sstring_Res3 *res, long N, long n, int r, int s)
{

    sstring_Run(gen, res, N, n, r, s);
    std::string pVal = std::to_string(res->NRuns->pVal2[gofw_Mean]) + ",";
    pVal += std::to_string(res->NBits->pVal2[gofw_Mean]);
    return pVal;
}

std::string sstring_AutoCor_nu(unif01_Gen *gen, sres_Basic *res, long N, long n, int r, int s, int d)
{

    sstring_AutoCor(gen, res, 10, 30 + BILLION, r, s, 1);
    std::string pVal = std::to_string(res->pVal2[gofw_Sum]);
    return pVal;
}

std::string GetPVal_Walk_nu(long N, swalk_Res *res)
{
    if (N == 1)
    {
        std::string pVal = std::to_string(res->H[0]->pVal2[gofw_Mean]) + ",";
        pVal += std::to_string(res->M[0]->pVal2[gofw_Mean]) + ",";
        pVal += std::to_string(res->J[0]->pVal2[gofw_Mean]) + ",";
        pVal += std::to_string(res->R[0]->pVal2[gofw_Mean]) + ",";
        pVal += std::to_string(res->C[0]->pVal2[gofw_Mean]);

        return pVal;
    }
    else
    {
        std::string pVal = std::to_string(res->H[0]->pVal2[gofw_Sum]) + ",";
        pVal += std::to_string(res->M[0]->pVal2[gofw_Sum]) + ",";
        pVal += std::to_string(res->J[0]->pVal2[gofw_Sum]) + ",";
        pVal += std::to_string(res->R[0]->pVal2[gofw_Sum]) + ",";
        pVal += std::to_string(res->C[0]->pVal2[gofw_Sum]);

        return pVal;
    }
}

std::string GetPVal_CPairs_nu(long N, snpair_Res *res)
/* Get the p-values in a snpair_ClosePairs test */
{
    if (N == 1)
    {
        std::string pVal = std::to_string(res->pVal[snpair_NP]) + ",";

        pVal += std::to_string(res->pVal[snpair_mNP]) + "," + "," + "," + ",";
        return pVal;
    }
    else
    {
        std::string pVal = std::to_string(res->pVal[snpair_NP]) + ",";

        pVal += std::to_string(res->pVal[snpair_mNP]) + ",";

        pVal += std::to_string(res->pVal[snpair_mNP1]) + ",";

        pVal += std::to_string(res->pVal[snpair_mNP2]) + ",";

        pVal += std::to_string(res->pVal[snpair_NJumps]);

        if (snpair_mNP2S_Flag)
        {
            pVal += "," + std::to_string(res->pVal[snpair_mNP2S]);
        }
        else
        {
            pVal += ",";
        }
        return pVal;
    }
}

std::string selectSmallCrushTest(unif01_Gen *gen, int testNum)
{
    swrite_Basic = FALSE;
    swrite_Host = FALSE;
    const int r = 0;

    switch (testNum)
    {
    case 1:
    {
        sres_Poisson *res1 = sres_CreatePoisson();
        std::string pVal = smarsa_BirthdaySpacings_nu(gen, res1, 1, 5 * MILLION, r, 1073741824, 2, 1);
        sres_DeletePoisson(res1);

        return pVal;
    }
    case 2:
    {
        sknuth_Res2 *res2 = sknuth_CreateRes2();
        std::string pVal = sknuth_Collision_nu(gen, res2, 1, 5 * MILLION, 0, 65536, 2);
        sknuth_DeleteRes2(res2);

        return pVal;
    }
    case 3:
    {
        sres_Chi2 *res3 = sres_CreateChi2();
        std::string pVal = sknuth_Gap_nu(gen, res3, 1, MILLION / 5, 22, 0.0, .00390625);
        sres_DeleteChi2(res3);

        return pVal;
    }
    case 4:
    {
        sres_Chi2 *res4 = sres_CreateChi2();
        std::string pVal = sknuth_SimpPoker_nu(gen, res4, 1, 2 * MILLION / 5, 24, 64, 64);
        sres_DeleteChi2(res4);

        return pVal;
    }
    case 5:
    {
        sres_Chi2 *res5 = sres_CreateChi2();
        std::string pVal = sknuth_CouponCollector_nu(gen, res5, 1, MILLION / 2, 26, 16);
        sres_DeleteChi2(res5);

        return pVal;
    }
    case 6:
    {
        sknuth_Res1 *res6 = sknuth_CreateRes1();
        std::string pVal = sknuth_MaxOft_nu_Mean_Mean(gen, res6, 1, 2 * MILLION, 0, MILLION / 10, 6);
        sknuth_DeleteRes1(res6);

        return pVal;
    }
    case 7:
    {
        sres_Chi2 *res7 = sres_CreateChi2();
        std::string pVal = svaria_WeightDistrib_nu(gen, res7, 1, MILLION / 5, 27, 256, 0.0, 0.125);
        sres_DeleteChi2(res7);

        return pVal;
    }
    case 8:
    {
        sres_Chi2 *res8 = sres_CreateChi2();
        std::string pVal = smarsa_MatrixRank_nu_Mean(gen, res8, 1, 20 * THOUSAND, 20, 10, 60, 60);
        sres_DeleteChi2(res8);

        return pVal;
    }
    case 9:
    {
        sstring_Res *res9 = sstring_CreateRes();
        std::string pVal = sstring_HammingIndep_nu_Mean(gen, res9, 1, MILLION / 2, 20, 10, 300, 0);
        sstring_DeleteRes(res9);

        return pVal;
    }
    case 10:
    {
        swalk_Res *res10 = swalk_CreateRes();
        std::string pVal = swalk_RandomWalk1_nu(gen, res10, 1, MILLION, r, 30, 150, 150, 1);
        swalk_DeleteRes(res10);

        return pVal;
    }
    default:
    {
        std::string pVal = "Invalid case.\n";

        return pVal;
    }
    }
}

std::string selectBigCrushTest(unif01_Gen *gen, int testNum)
{
    swrite_Basic = FALSE;
    swrite_Host = FALSE;
    const int s = 30;
    const int r = 0;
    lebool flag = snpair_mNP2S_Flag;
    switch (testNum)
    {
    case 1:
    {
        sres_Basic *res1 = sres_CreateBasic();
        std::string pVal = smarsa_SerialOver_nu(gen, res1, 1, BILLION, 0, 256, 3);
        sres_DeleteBasic(res1);

        return pVal;
    }
    case 2:
    {
        sres_Basic *res2 = sres_CreateBasic();
        std::string pVal = smarsa_SerialOver_nu(gen, res2, 1, BILLION, 22, 256, 3);
        sres_DeleteBasic(res2);

        return pVal;
    }
    case 3:
    {
        smarsa_Res *res3 = smarsa_CreateRes();
        std::string pVal = smarsa_CollisionOver_nu(gen, res3, 30, 20 * MILLION, 0, 1024 * 1024 * 2, 2);
        smarsa_DeleteRes(res3);

        return pVal;
    }
    case 4:
    {
        smarsa_Res *res4 = smarsa_CreateRes();
        std::string pVal = smarsa_CollisionOver_nu(gen, res4, 30, 20 * MILLION, 9, 1024 * 1024 * 2, 2);
        smarsa_DeleteRes(res4);

        return pVal;
    }
    case 5:
    {
        smarsa_Res *res5 = smarsa_CreateRes();
        std::string pVal = smarsa_CollisionOver_nu(gen, res5, 30, 20 * MILLION, 0, 1024 * 16, 3);
        smarsa_DeleteRes(res5);

        return pVal;
    }
    case 6:
    {
        smarsa_Res *res6 = smarsa_CreateRes();
        std::string pVal = smarsa_CollisionOver_nu(gen, res6, 30, 20 * MILLION, 16, 1024 * 16, 3);
        smarsa_DeleteRes(res6);

        return pVal;
    }
    case 7:
    {
        smarsa_Res *res7 = smarsa_CreateRes();
        std::string pVal = smarsa_CollisionOver_nu(gen, res7, 30, 20 * MILLION, 0, 64, 7);
        smarsa_DeleteRes(res7);

        return pVal;
    }
    case 8:
    {
        smarsa_Res *res8 = smarsa_CreateRes();
        std::string pVal = smarsa_CollisionOver_nu(gen, res8, 30, 20 * MILLION, 24, 64, 7);
        smarsa_DeleteRes(res8);

        return pVal;
    }
    case 9:
    {
        smarsa_Res *res9 = smarsa_CreateRes();
        std::string pVal = smarsa_CollisionOver_nu(gen, res9, 30, 20 * MILLION, 0, 8, 14);
        smarsa_DeleteRes(res9);

        return pVal;
    }
    case 10:
    {
        smarsa_Res *res10 = smarsa_CreateRes();
        std::string pVal = smarsa_CollisionOver_nu(gen, res10, 30, 20 * MILLION, 27, 8, 14);
        smarsa_DeleteRes(res10);

        return pVal;
    }
    case 11:
    {
        smarsa_Res *res11 = smarsa_CreateRes();
        std::string pVal = smarsa_CollisionOver_nu(gen, res11, 30, 20 * MILLION, 0, 4, 21);
        smarsa_DeleteRes(res11);

        return pVal;
    }
    case 12:
    {
        smarsa_Res *res12 = smarsa_CreateRes();
        std::string pVal = smarsa_CollisionOver_nu(gen, res12, 30, 20 * MILLION, 28, 4, 21);
        smarsa_DeleteRes(res12);

        return pVal;
    }
    case 13:
    {
        sres_Poisson *res13 = sres_CreatePoisson();
        std::string pVal = "";
#ifdef USE_LONGLONG
        long d;
#if LONG_MAX <= 2147483647L
        d = 1073741824L;
        pVal += smarsa_BirthdaySpacings_nu(gen, res13, 250, 4 * MILLION, 0, d, 2, 1);
#else
        d = 2147483648L;
        pVal += smarsa_BirthdaySpacings_nu(gen, res13, 100, 10 * MILLION, 0, d, 2, 1);
#endif
#else
        pVal += smnarsa_BirthdaySpacings_nu(gen, res13, 10 * THOUSAND, MILLION / 10, 0, 67108864, 2, 1);
#endif
        sres_DeletePoisson(res13);

        return pVal;
    }
    case 14:
    {
        sres_Poisson *res14 = sres_CreatePoisson();
        std::string pVal = "";
#ifdef USE_LONGLONG
        pVal += smarsa_BirthdaySpacings_nu(gen, res14, 20, 20 * MILLION, 0, 2097152, 3, 1);
#else
        pVal += smarsa_BirthdaySpacings_nu(gen, res14, 10 * THOUSAND, MILLION / 10, 0, 1024 * 8, 4, 1);
#endif
        sres_DeletePoisson(res14);

        return pVal;
    }
    case 15:
    {
        sres_Poisson *res15 = sres_CreatePoisson();
        std::string pVal = "";
#ifdef USE_LONGLONG
        pVal += smarsa_BirthdaySpacings_nu(gen, res15, 20, 30 * MILLION, 14, 65536, 4, 1);
#else
        pVal += smarsa_BirthdaySpacings_nu(gen, res15, 10 * THOUSAND, MILLION / 10, 16, 1024 * 8, 4, 1);
#endif
        sres_DeletePoisson(res15);

        return pVal;
    }
    case 16:
    {
        sres_Poisson *res16 = sres_CreatePoisson();
        std::string pVal = "";
#ifdef USE_LONGLONG
        pVal += smarsa_BirthdaySpacings_nu(gen, res16, 20, 20 * MILLION, 0, 512, 7, 1);
#else
        pVal += smarsa_BirthdaySpacings_nu(gen, res16, 10 * THOUSAND, MILLION / 10, 0, 16, 13, 1);
#endif
        sres_DeletePoisson(res16);

        return pVal;
    }
    case 17:
    {
        sres_Poisson *res17 = sres_CreatePoisson();
        std::string pVal = "";
#ifdef USE_LONGLONG
        pVal += smarsa_BirthdaySpacings_nu(gen, res17, 20, 20 * MILLION, 7, 512, 7, 1);
#else
        pVal += smarsa_BirthdaySpacings_nu(gen, res17, 10 * THOUSAND, MILLION / 10, 5, 16, 13, 1);
#endif
        sres_DeletePoisson(res17);

        return pVal;
    }
    case 18:
    {
        sres_Poisson *res18 = sres_CreatePoisson();
        std::string pVal = "";
#ifdef USE_LONGLONG
        pVal += smarsa_BirthdaySpacings_nu(gen, res18, 20, 30 * MILLION, 14, 256, 8, 1);
#else
        pVal += smarsa_BirthdaySpacings_nu(gen, res18, 10 * THOUSAND, MILLION / 10, 10, 16, 13, 1);
#endif
        sres_DeletePoisson(res18);

        return pVal;
    }
    case 19:
    {
        sres_Poisson *res19 = sres_CreatePoisson();
        std::string pVal = "";
#ifdef USE_LONGLONG
        pVal += smarsa_BirthdaySpacings_nu(gen, res19, 20, 30 * MILLION, 22, 256, 8, 1);
#else
        pVal += smarsa_BirthdaySpacings_nu(gen, res19, 10 * THOUSAND, MILLION / 10, 15, 16, 13, 1);
#endif
        sres_DeletePoisson(res19);

        return pVal;
    }
    case 20:
    {
        sres_Poisson *res20 = sres_CreatePoisson();
        std::string pVal = "";
#ifdef USE_LONGLONG
        pVal += smarsa_BirthdaySpacings_nu(gen, res20, 20, 30 * MILLION, 0, 16, 16, 1);
#else
        pVal += smarsa_BirthdaySpacings_nu(gen, res20, 10 * THOUSAND, MILLION / 10, 20, 16, 13, 1);
#endif
        sres_DeletePoisson(res20);

        return pVal;
    }
    case 21:
    {
        sres_Poisson *res21 = sres_CreatePoisson();
        std::string pVal = "";
#ifdef USE_LONGLONG
        pVal += smarsa_BirthdaySpacings_nu(gen, res21, 20, 30 * MILLION, 26, 16, 16, 1);
#else
        pVal += smarsa_BirthdaySpacings_nu(gen, res21, 10 * THOUSAND, MILLION / 10, 26, 16, 13, 1);
#endif
        sres_DeletePoisson(res21);

        return pVal;
    }
    case 22:
    {
        snpair_Res *res22 = snpair_CreateRes();
        lebool flag = snpair_mNP2S_Flag;
        snpair_mNP2S_Flag = TRUE;
        std::string pVal = snpair_ClosePairs_GetPVal_CPairs_nu(gen, res22, 30, 6 * MILLION, 0, 3, 0, 30, 40);
        snpair_DeleteRes(res22);
        snpair_mNP2S_Flag = flag;

        return pVal;
    }
    case 23:
    {
        snpair_Res *res23 = snpair_CreateRes();
        lebool flag = snpair_mNP2S_Flag;
        snpair_mNP2S_Flag = TRUE;
        std::string pVal = snpair_ClosePairs_GetPVal_CPairs_nu(gen, res23, 20, 4 * MILLION, 0, 5, 0, 30, 40);
        snpair_DeleteRes(res23);
        snpair_mNP2S_Flag = flag;

        return pVal;
    }
    case 24:
    {
        snpair_Res *res24 = snpair_CreateRes();
        lebool flag = snpair_mNP2S_Flag;
        snpair_mNP2S_Flag = TRUE;
        std::string pVal = snpair_ClosePairs_GetPVal_CPairs_nu(gen, res24, 10, 3 * MILLION, 0, 9, 0, 30, 20);
        snpair_DeleteRes(res24);
        snpair_mNP2S_Flag = flag;

        return pVal;
    }
    case 25:
    {
        snpair_Res *res25 = snpair_CreateRes();
        lebool flag = snpair_mNP2S_Flag;
        snpair_mNP2S_Flag = TRUE;
        std::string pVal = snpair_ClosePairs_GetPVal_CPairs_nu(gen, res25, 5, 2 * MILLION, 0, 16, 0, 30, 10);
        snpair_DeleteRes(res25);
        snpair_mNP2S_Flag = flag;

        return pVal;
    }
    case 26:
    {
        sres_Chi2 *res26 = sres_CreateChi2();
        std::string pVal = sknuth_SimpPoker_nu(gen, res26, 1, 400 * MILLION, 0, 8, 8);
        sres_DeleteChi2(res26);

        return pVal;
    }
    case 27:
    {
        sres_Chi2 *res27 = sres_CreateChi2();
        std::string pVal = sknuth_SimpPoker_nu(gen, res27, 1, 400 * MILLION, 27, 8, 8);
        sres_DeleteChi2(res27);

        return pVal;
    }
    case 28:
    {
        sres_Chi2 *res28 = sres_CreateChi2();
        std::string pVal = sknuth_SimpPoker_nu(gen, res28, 1, 100 * MILLION, 0, 32, 32);
        sres_DeleteChi2(res28);

        return pVal;
    }
    case 29:
    {
        sres_Chi2 *res29 = sres_CreateChi2();
        std::string pVal = sknuth_SimpPoker_nu(gen, res29, 1, 100 * MILLION, 25, 32, 32);
        sres_DeleteChi2(res29);

        return pVal;
    }
    case 30:
    {
        sres_Chi2 *res30 = sres_CreateChi2();
        std::string pVal = sknuth_CouponCollector_nu(gen, res30, 1, 200 * MILLION, 0, 8);
        sres_DeleteChi2(res30);

        return pVal;
    }
    case 31:
    {
        sres_Chi2 *res31 = sres_CreateChi2();
        std::string pVal = sknuth_CouponCollector_nu(gen, res31, 1, 200 * MILLION, 10, 8);
        sres_DeleteChi2(res31);

        return pVal;
    }
    case 32:
    {
        sres_Chi2 *res32 = sres_CreateChi2();
        std::string pVal = sknuth_CouponCollector_nu(gen, res32, 1, 200 * MILLION, 20, 8);
        sres_DeleteChi2(res32);

        return pVal;
    }
    case 33:
    {
        sres_Chi2 *res33 = sres_CreateChi2();
        std::string pVal = sknuth_CouponCollector_nu(gen, res33, 1, 200 * MILLION, 27, 8);
        sres_DeleteChi2(res33);

        return pVal;
    }
    case 34:
    {
        sres_Chi2 *res34 = sres_CreateChi2();
        std::string pVal = sknuth_Gap_nu(gen, res34, 1, BILLION / 2, 0, 0.0, 1.0 / 16.0);
        sres_DeleteChi2(res34);

        return pVal;
    }
    case 35:
    {
        sres_Chi2 *res35 = sres_CreateChi2();
        std::string pVal = sknuth_Gap_nu(gen, res35, 1, 300 * MILLION, 25, 0.0, 1.0 / 32.0);
        sres_DeleteChi2(res35);

        return pVal;
    }
    case 36:
    {
        sres_Chi2 *res36 = sres_CreateChi2();
        std::string pVal = sknuth_Gap_nu(gen, res36, 1, BILLION / 10, 0, 0.0, 1.0 / 128.0);
        sres_DeleteChi2(res36);

        return pVal;
    }
    case 37:
    {
        sres_Chi2 *res37 = sres_CreateChi2();
        std::string pVal = sknuth_Gap_nu(gen, res37, 1, 10 * MILLION, 20, 0.0, 1.0 / 1024.0);
        sres_DeleteChi2(res37);

        return pVal;
    }
    case 38:
    {
        sres_Chi2 *res38 = sres_CreateChi2();
        std::string pVal = sknuth_Run_nu_Sum(gen, res38, 5, BILLION, 0, FALSE);
        sres_DeleteChi2(res38);

        return pVal;
    }
    case 39:
    {
        sres_Chi2 *res39 = sres_CreateChi2();
        std::string pVal = sknuth_Run_nu_Sum(gen, res39, 10, BILLION, 15, TRUE);
        sres_DeleteChi2(res39);

        return pVal;
    }
    case 40:
    {
        sres_Chi2 *res40 = sres_CreateChi2();
        std::string pVal = sknuth_Permutation_nu(gen, res40, 1, BILLION, 5, 3);
        sres_DeleteChi2(res40);

        return pVal;
    }
    case 41:
    {
        sres_Chi2 *res41 = sres_CreateChi2();
        std::string pVal = sknuth_Permutation_nu(gen, res41, 1, BILLION, 5, 5);
        sres_DeleteChi2(res41);

        return pVal;
    }
    case 42:
    {
        sres_Chi2 *res42 = sres_CreateChi2();
        std::string pVal = sknuth_Permutation_nu(gen, res42, 1, BILLION / 2, 5, 7);
        sres_DeleteChi2(res42);

        return pVal;
    }
    case 43:
    {
        sres_Chi2 *res43 = sres_CreateChi2();
        std::string pVal = sknuth_Permutation_nu(gen, res43, 1, BILLION / 2, 10, 10);
        sres_DeleteChi2(res43);

        return pVal;
    }
    case 44:
    {
        sknuth_Res2 *res44 = sknuth_CreateRes2();
        std::string pVal = sknuth_CollisionPermut_nu(gen, res44, 20, 20 * MILLION, 0, 14);
        sknuth_DeleteRes2(res44);

        return pVal;
    }
    case 45:
    {
        sknuth_Res2 *res45 = sknuth_CreateRes2();
        std::string pVal = sknuth_CollisionPermut_nu(gen, res45, 20, 20 * MILLION, 10, 14);
        sknuth_DeleteRes2(res45);

        return pVal;
    }
    case 46:
    {
        sknuth_Res1 *res46 = sknuth_CreateRes1();
        std::string pVal = sknuth_MaxOft_nu_Sum_AD(gen, res46, 40, 10 * MILLION, 0, MILLION / 10, 8);
        sknuth_DeleteRes1(res46);

        return pVal;
    }
    case 47:
    {
        sknuth_Res1 *res47 = sknuth_CreateRes1();
        std::string pVal = sknuth_MaxOft_nu_Sum_AD(gen, res47, 30, 10 * MILLION, 0, MILLION / 10, 16);
        sknuth_DeleteRes1(res47);

        return pVal;
    }
    case 48:
    {
        sknuth_Res1 *res48 = sknuth_CreateRes1();
        std::string pVal = sknuth_MaxOft_nu_Sum_AD(gen, res48, 20, 10 * MILLION, 0, MILLION / 10, 24);
        sknuth_DeleteRes1(res48);

        return pVal;
    }
    case 49:
    {
        sknuth_Res1 *res49 = sknuth_CreateRes1();
        std::string pVal = sknuth_MaxOft_nu_Sum_AD(gen, res49, 20, 10 * MILLION, 0, MILLION / 10, 32);
        sknuth_DeleteRes1(res49);

        return pVal;
    }
    case 50:
    {
        sres_Basic *res50 = sres_CreateBasic();
        std::string pVal = svaria_SampleProd_nu_AD(gen, res50, 40, 10 * MILLION, 0, 8);
        sres_DeleteBasic(res50);

        return pVal;
    }
    case 51:
    {
        sres_Basic *res51 = sres_CreateBasic();
        std::string pVal = svaria_SampleProd_nu_AD(gen, res51, 20, 10 * MILLION, 0, 16);
        sres_DeleteBasic(res51);

        return pVal;
    }
    case 52:
    {
        sres_Basic *res52 = sres_CreateBasic();
        std::string pVal = svaria_SampleProd_nu_AD(gen, res52, 20, 10 * MILLION, 0, 24);
        sres_DeleteBasic(res52);

        return pVal;
    }
    case 53:
    {
        sres_Basic *res53 = sres_CreateBasic();
        std::string pVal = svaria_SampleMean_nu(gen, res53, 20 * MILLION, 30, 0);
        sres_DeleteBasic(res53);

        return pVal;
    }
    case 54:
    {
        sres_Basic *res54 = sres_CreateBasic();
        std::string pVal = svaria_SampleMean_nu(gen, res54, 20 * MILLION, 30, 10);
        sres_DeleteBasic(res54);

        return pVal;
    }
    case 55:
    {
        sres_Basic *res55 = sres_CreateBasic();
        std::string pVal = svaria_SampleCorr_nu(gen, res55, 1, 2 * BILLION, 0, 1);
        sres_DeleteBasic(res55);

        return pVal;
    }
    case 56:
    {
        sres_Basic *res56 = sres_CreateBasic();
        std::string pVal = svaria_SampleCorr_nu(gen, res56, 1, 2 * BILLION, 0, 2);
        sres_DeleteBasic(res56);

        return pVal;
    }
    case 57:
    {
        sres_Basic *res57 = sres_CreateBasic();
        std::string pVal = svaria_AppearanceSpacings_nu(gen, res57, 1, 10 * MILLION, BILLION, r, 3, 15);
        sres_DeleteBasic(res57);

        return pVal;
    }
    case 58:
    {
        sres_Basic *res58 = sres_CreateBasic();
        std::string pVal = svaria_AppearanceSpacings_nu(gen, res58, 1, 10 * MILLION, BILLION, 27, 3, 15);
        sres_DeleteBasic(res58);

        return pVal;
    }
    case 59:
    {
        sres_Chi2 *res59 = sres_CreateChi2();
        std::string pVal = svaria_WeightDistrib_nu(gen, res59, 1, 20 * MILLION, 0, 256, 0.0, 0.25);
        sres_DeleteChi2(res59);

        return pVal;
    }
    case 60:
    {
        sres_Chi2 *res60 = sres_CreateChi2();
        std::string pVal = svaria_WeightDistrib_nu(gen, res60, 1, 20 * MILLION, 20, 256, 0.0, 0.25);
        sres_DeleteChi2(res60);

        return pVal;
    }
    case 61:
    {
        sres_Chi2 *res61 = sres_CreateChi2();
        std::string pVal = svaria_WeightDistrib_nu(gen, res61, 1, 20 * MILLION, 28, 256, 0.0, 0.25);
        sres_DeleteChi2(res61);

        return pVal;
    }
    case 62:
    {
        sres_Chi2 *res62 = sres_CreateChi2();
        std::string pVal = svaria_WeightDistrib_nu(gen, res62, 1, 20 * MILLION, 0, 256, 0.0, 0.0625);
        sres_DeleteChi2(res62);

        return pVal;
    }
    case 63:
    {
        sres_Chi2 *res63 = sres_CreateChi2();
        std::string pVal = svaria_WeightDistrib_nu(gen, res63, 1, 20 * MILLION, 10, 256, 0.0, 0.0625);
        sres_DeleteChi2(res63);

        return pVal;
    }
    case 64:
    {
        sres_Chi2 *res64 = sres_CreateChi2();
        std::string pVal = svaria_WeightDistrib_nu(gen, res64, 1, 20 * MILLION, 26, 256, 0.0, 0.0625);
        sres_DeleteChi2(res64);

        return pVal;
    }
    case 65:
    {
        sres_Chi2 *res65 = sres_CreateChi2();
        std::string pVal = svaria_SumCollector_nu(gen, res65, 1, 500 * MILLION, 0, 10.0);
        sres_DeleteChi2(res65);

        return pVal;
    }
    case 66:
    {
        sres_Chi2 *res66 = sres_CreateChi2();
        std::string pVal = smarsa_MatrixRank_nu_Sum(gen, res66, 10, MILLION, r, 5, 30, 30);
        sres_DeleteChi2(res66);

        return pVal;
    }
    case 67:
    {
        sres_Chi2 *res67 = sres_CreateChi2();
        std::string pVal = smarsa_MatrixRank_nu_Sum(gen, res67, 10, MILLION, 25, 5, 30, 30);
        sres_DeleteChi2(res67);

        return pVal;
    }
    case 68:
    {
        sres_Chi2 *res68 = sres_CreateChi2();
        std::string pVal = smarsa_MatrixRank_nu_Mean(gen, res68, 1, 5 * THOUSAND, r, 4, 1000, 1000);
        sres_DeleteChi2(res68);

        return pVal;
    }
    case 69:
    {
        sres_Chi2 *res69 = sres_CreateChi2();
        std::string pVal = smarsa_MatrixRank_nu_Mean(gen, res69, 1, 5 * THOUSAND, 26, 4, 1000, 1000);
        sres_DeleteChi2(res69);

        return pVal;
    }
    case 70:
    {
        sres_Chi2 *res70 = sres_CreateChi2();
        std::string pVal = smarsa_MatrixRank_nu_Mean(gen, res70, 1, 80, 15, 15, 5000, 5000);
        sres_DeleteChi2(res70);

        return pVal;
    }
    case 71:
    {
        sres_Chi2 *res71 = sres_CreateChi2();
        std::string pVal = smarsa_MatrixRank_nu_Mean(gen, res71, 1, 80, 0, 30, 5000, 5000);
        sres_DeleteChi2(res71);

        return pVal;
    }
    case 72:
    {
        sres_Chi2 *res72 = sres_CreateChi2();
        std::string pVal = smarsa_Savir2_nu_Sum(gen, res72, 10, 10 * MILLION, 10, 1024 * 1024, 30);
        sres_DeleteChi2(res72);

        return pVal;
    }
    case 73:
    {
        smarsa_Res2 *res73 = smarsa_CreateRes2();
        std::string pVal = smarsa_GCD_nu_Sum(gen, res73, 10, 50 * MILLION, 0, 30);
        smarsa_DeleteRes2(res73);

        return pVal;
    }
    case 74:
    {
        swalk_Res *res74 = swalk_CreateRes();
        std::string pVal = swalk_RandomWalk1_nu(gen, res74, 1, 100 * MILLION, r, 5, 50, 50, 1);
        swalk_DeleteRes(res74);

        return pVal;
    }
    case 75:
    {
        swalk_Res *res75 = swalk_CreateRes();
        std::string pVal = swalk_RandomWalk1_nu(gen, res75, 1, 100 * MILLION, 25, 5, 50, 50, 1);
        swalk_DeleteRes(res75);

        return pVal;
    }
    case 76:
    {
        swalk_Res *res76 = swalk_CreateRes();
        std::string pVal = swalk_RandomWalk1_nu(gen, res76, 1, 10 * MILLION, r, 10, 1000, 1000, 1);
        swalk_DeleteRes(res76);

        return pVal;
    }
    case 77:
    {
        swalk_Res *res77 = swalk_CreateRes();
        std::string pVal = swalk_RandomWalk1_nu(gen, res77, 1, 10 * MILLION, 20, 10, 1000, 1000, 1);
        swalk_DeleteRes(res77);

        return pVal;
    }
    case 78:
    {
        swalk_Res *res78 = swalk_CreateRes();
        std::string pVal = swalk_RandomWalk1_nu(gen, res78, 1, 1 * MILLION, r, 15, 10000, 10000, 1);
        swalk_DeleteRes(res78);

        return pVal;
    }
    case 79:
    {
        swalk_Res *res79 = swalk_CreateRes();
        std::string pVal = swalk_RandomWalk1_nu(gen, res79, 1, 1 * MILLION, 15, 15, 10000, 10000, 1);
        swalk_DeleteRes(res79);

        return pVal;
    }
    case 80:
    {
        scomp_Res *res80 = scomp_CreateRes();
        std::string pVal = scomp_LinearComp_nu(gen, res80, 1, 400 * THOUSAND + 20, r, 1);
        scomp_DeleteRes(res80);

        return pVal;
    }
    case 81:
    {
        scomp_Res *res81 = scomp_CreateRes();
        std::string pVal = scomp_LinearComp_nu(gen, res81, 1, 400 * THOUSAND + 20, 29, 1);
        scomp_DeleteRes(res81);

        return pVal;
    }
    case 82:
    {
        sres_Basic *res82 = sres_CreateBasic();
        std::string pVal = scomp_LempelZiv_nu(gen, res82, 10, 27, r, s);
        sres_DeleteBasic(res82);

        return pVal;
    }
    case 83:
    {
        sres_Basic *res83 = sres_CreateBasic();
        std::string pVal = scomp_LempelZiv_nu(gen, res83, 10, 27, 15, 15);
        sres_DeleteBasic(res83);

        return pVal;
    }
    case 84:
    {
        sspectral_Res *res84 = sspectral_CreateRes();
        std::string pVal = sspectral_Fourier3_nu(gen, res84, 100 * THOUSAND, 14, r, 3);
        sspectral_DeleteRes(res84);

        return pVal;
    }
    case 85:
    {
        sspectral_Res *res85 = sspectral_CreateRes();
        std::string pVal = sspectral_Fourier3_nu(gen, res85, 100 * THOUSAND, 14, 27, 3);
        sspectral_DeleteRes(res85);

        return pVal;
    }
    case 86:
    {
        sstring_Res2 *res86 = sstring_CreateRes2();
        std::string pVal = sstring_LongestHeadRun_nu(gen, res86, 1, 1000, r, 3, 20 + 10 * MILLION);
        sstring_DeleteRes2(res86);

        return pVal;
    }
    case 87:
    {
        sstring_Res2 *res87 = sstring_CreateRes2();
        std::string pVal = sstring_LongestHeadRun_nu(gen, res87, 1, 1000, 27, 3, 20 + 10 * MILLION);
        sstring_DeleteRes2(res87);

        return pVal;
    }
    case 88:
    {
        sres_Chi2 *res88 = sres_CreateChi2();
        std::string pVal = sstring_PeriodsInStrings_nu_Sum(gen, res88, 10, BILLION / 2, r, 10);
        sres_DeleteChi2(res88);

        return pVal;
    }
    case 89:
    {
        sres_Chi2 *res89 = sres_CreateChi2();
        std::string pVal = sstring_PeriodsInStrings_nu_Sum(gen, res89, 10, BILLION / 2, r, 10);
        sres_DeleteChi2(res89);

        return pVal;
    }
    case 90:
    {
        sres_Basic *res90 = sres_CreateBasic();
        std::string pVal = sstring_HammingWeight2_nu(gen, res90, 10, BILLION, r, 3, MILLION);
        sres_DeleteBasic(res90);

        return pVal;
    }
    case 91:
    {
        sres_Basic *res91 = sres_CreateBasic();
        std::string pVal = sstring_HammingWeight2_nu(gen, res91, 10, BILLION, 27, 3, MILLION);
        sres_DeleteBasic(res91);

        return pVal;
    }
    case 92:
    {
        sstring_Res *res92 = sstring_CreateRes();
        std::string pVal = sstring_HammingCorr_nu(gen, res92, 1, BILLION, 10, 10, s);
        sstring_DeleteRes(res92);

        return pVal;
    }
    case 93:
    {
        sstring_Res *res93 = sstring_CreateRes();
        std::string pVal = sstring_HammingCorr_nu(gen, res93, 1, 100 * MILLION, 10, 10, 10 * s);
        sstring_DeleteRes(res93);

        return pVal;
    }
    case 94:
    {
        sstring_Res *res94 = sstring_CreateRes();
        std::string pVal = sstring_HammingCorr_nu(gen, res94, 1, 100 * MILLION, 10, 10, 40 * s);
        sstring_DeleteRes(res94);

        return pVal;
    }
    case 95:
    {
        sstring_Res *res95 = sstring_CreateRes();
        std::string pVal = sstring_HammingIndep_nu_Sum(gen, res95, 10, 30 * MILLION, r, 3, s, 0);
        sstring_DeleteRes(res95);

        return pVal;
    }
    case 96:
    {
        sstring_Res *res96 = sstring_CreateRes();
        std::string pVal = sstring_HammingIndep_nu_Sum(gen, res96, 10, 30 * MILLION, 27, 3, s, 0);
        sstring_DeleteRes(res96);

        return pVal;
    }
    case 97:
    {
        sstring_Res *res97 = sstring_CreateRes();
        std::string pVal = sstring_HammingIndep_nu_Mean(gen, res97, 1, 30 * MILLION, r, 4, 10 * s, 0);
        sstring_DeleteRes(res97);

        return pVal;
    }
    case 98:
    {
        sstring_Res *res98 = sstring_CreateRes();
        std::string pVal = sstring_HammingIndep_nu_Mean(gen, res98, 1, 30 * MILLION, 26, 4, 10 * s, 0);
        sstring_DeleteRes(res98);

        return pVal;
    }
    case 99:
    {
        sstring_Res *res99 = sstring_CreateRes();
        std::string pVal = sstring_HammingIndep_nu_Mean(gen, res99, 1, 10 * MILLION, r, 5, 40 * s, 0);
        sstring_DeleteRes(res99);

        return pVal;
    }
    case 100:
    {
        sstring_Res *res100 = sstring_CreateRes();
        std::string pVal = sstring_HammingIndep_nu_Mean(gen, res100, 1, 10 * MILLION, 25, 5, 40 * s, 0);
        sstring_DeleteRes(res100);

        return pVal;
    }
    case 101:
    {
        sstring_Res3 *res101 = sstring_CreateRes3();
        std::string pVal = sstring_Run_nu(gen, res101, 1, 2 * BILLION, r, 3);
        sstring_DeleteRes3(res101);

        return pVal;
    }
    case 102:
    {
        sstring_Res3 *res102 = sstring_CreateRes3();
        std::string pVal = sstring_Run_nu(gen, res102, 1, 2 * BILLION, 27, 3);
        sstring_DeleteRes3(res102);

        return pVal;
    }
    case 103:
    {
        sres_Basic *res103 = sres_CreateBasic();
        std::string pVal = sstring_AutoCor_nu(gen, res103, 10, 30 + BILLION, r, 3, 1);
        sres_DeleteBasic(res103);

        return pVal;
    }
    case 104:
    {
        sres_Basic *res104 = sres_CreateBasic();
        std::string pVal = sstring_AutoCor_nu(gen, res104, 10, 30 + BILLION, r, 3, 3);
        sres_DeleteBasic(res104);

        return pVal;
    }
    case 105:
    {
        sres_Basic *res105 = sres_CreateBasic();
        std::string pVal = sstring_AutoCor_nu(gen, res105, 10, 30 + BILLION, 27, 3, 1);
        sres_DeleteBasic(res105);

        return pVal;
    }
    case 106:
    {
        sres_Basic *res106 = sres_CreateBasic();
        std::string pVal = sstring_AutoCor_nu(gen, res106, 10, 30 + BILLION, 27, 3, 3);
        sres_DeleteBasic(res106);

        return pVal;
    }
    default:
    {
        std::string pVal = "Invalid case\n";
        return pVal;
    }
    }
}

std::string selectCrushTest(unif01_Gen *gen, int testNum)
{
    swrite_Basic = FALSE;
    swrite_Host = FALSE;
    const int s = 30;
    const int r = 0;
    int i;

    switch (testNum)
    {
    case 1:
    {
        sres_Basic *res1 = sres_CreateBasic();
        std::string pVal = smarsa_SerialOver_nu(gen, res1, 1, 500 * MILLION, 0, 4096, 2);
        sres_DeleteBasic(res1);

        return pVal;
    }
    case 2:
    {
        sres_Basic *res2 = sres_CreateBasic();
        std::string pVal = smarsa_SerialOver_nu(gen, res2, 1, 300 * MILLION, 0, 64, 4);
        sres_DeleteBasic(res2);

        return pVal;
    }
    case 3:
    {
        smarsa_Res *res3 = smarsa_CreateRes();
        std::string pVal = smarsa_CollisionOver_nu(gen, res3, 10, 10 * MILLION, 0, 1024 * 1024, 2);
        smarsa_DeleteRes(res3);

        return pVal;
    }
    case 4:
    {
        smarsa_Res *res4 = smarsa_CreateRes();
        std::string pVal = smarsa_CollisionOver_nu(gen, res4, 10, 10 * MILLION, 10, 1024 * 1024, 2);
        smarsa_DeleteRes(res4);

        return pVal;
    }
    case 5:
    {
        smarsa_Res *res5 = smarsa_CreateRes();
        std::string pVal = smarsa_CollisionOver_nu(gen, res5, 10, 10 * MILLION, 0, 1024, 4);
        smarsa_DeleteRes(res5);

        return pVal;
    }
    case 6:
    {
        smarsa_Res *res6 = smarsa_CreateRes();
        std::string pVal = smarsa_CollisionOver_nu(gen, res6, 10, 10 * MILLION, 20, 1024, 4);
        smarsa_DeleteRes(res6);

        return pVal;
    }
    case 7:
    {
        smarsa_Res *res7 = smarsa_CreateRes();
        std::string pVal = smarsa_CollisionOver_nu(gen, res7, 10, 10 * MILLION, 0, 32, 8);
        smarsa_DeleteRes(res7);

        return pVal;
    }
    case 8:
    {
        smarsa_Res *res8 = smarsa_CreateRes();
        std::string pVal = smarsa_CollisionOver_nu(gen, res8, 10, 10 * MILLION, 25, 32, 8);
        smarsa_DeleteRes(res8);

        return pVal;
    }
    case 9:
    {
        smarsa_Res *res9 = smarsa_CreateRes();
        std::string pVal = smarsa_CollisionOver_nu(gen, res9, 10, 10 * MILLION, 0, 4, 20);
        smarsa_DeleteRes(res9);

        return pVal;
    }
    case 10:
    {
        smarsa_Res *res10 = smarsa_CreateRes();
        std::string pVal = smarsa_CollisionOver_nu(gen, res10, 10, 10 * MILLION, 28, 4, 20);
        smarsa_DeleteRes(res10);

        return pVal;
    }
    case 11:
    {
        sres_Poisson *res11 = sres_CreatePoisson();
        std::string pVal = "";
#ifdef USE_LONGLONG
        long d;
#if LONG_MAX <= 2147483647L
        d = 1073741824L;
        pVal += smarsa_BirthdaySpacings_nu(gen, res11, 10, 10 * MILLION, 0, d, 2, 1);
#else
        d = 2 * 1073741824L;
        pVal += smarsa_BirthdaySpacings_nu(gen, res11, 5, 20 * MILLION, 0, d, 2, 1);
#endif
#else
        pVal += smarsa_BirthdaySpacings_nu(gen, res11, 200, 4 * MILLION / 10, 0, 67108864, 2, 1);
#endif
        sres_DeletePoisson(res11);

        return pVal;
    }
    case 12:
    {
        sres_Poisson *res12 = sres_CreatePoisson();
        std::string pVal = "";
#ifdef USE_LONGLONG
        pVal += smarsa_BirthdaySpacings_nu(gen, res12, 5, 20 * MILLION, 0, 2097152, 3, 1);
#else
        pVal += smarsa_BirthdaySpacings_nu(gen, res12, 100, 4 * MILLION / 10, 0, 131072, 3, 1);
#endif
        sres_DeletePoisson(res12);

        return pVal;
    }
    case 13:
    {
        sres_Poisson *res13 = sres_CreatePoisson();
        std::string pVal = "";
#ifdef USE_LONGLONG
        pVal += smarsa_BirthdaySpacings_nu(gen, res13, 5, 20 * MILLION, 0, 65536, 4, 1);
#else
        pVal += smarsa_BirthdaySpacings_nu(gen, res13, 200, 4 * MILLION / 10, 0, 1024 * 8, 4, 1);
#endif
        sres_DeletePoisson(res13);

        return pVal;
    }
    case 14:
    {
        sres_Poisson *res14 = sres_CreatePoisson();
        std::string pVal = "";
#ifdef USE_LONGLONG
        pVal += smarsa_BirthdaySpacings_nu(gen, res14, 3, 20 * MILLION, 0, 512, 7, 1);
#else
        pVal += smarsa_BirthdaySpacings_nu(gen, res14, 100, 4 * MILLION / 10, 0, 16, 13, 1);
#endif
        sres_DeletePoisson(res14);

        return pVal;
    }
    case 15:
    {
        sres_Poisson *res15 = sres_CreatePoisson();
        std::string pVal = "";
#ifdef USE_LONGLONG
        pVal += smarsa_BirthdaySpacings_nu(gen, res15, 3, 20 * MILLION, 7, 512, 7, 1);
#else
        pVal += smarsa_BirthdaySpacings_nu(gen, res15, 100, 4 * MILLION / 10, 10, 16, 13, 1);
#endif
        sres_DeletePoisson(res15);

        return pVal;
    }
    case 16:
    {
        sres_Poisson *res16 = sres_CreatePoisson();
        std::string pVal = "";
#ifdef USE_LONGLONG
        pVal += smarsa_BirthdaySpacings_nu(gen, res16, 3, 20 * MILLION, 14, 256, 8, 1);
#else
        pVal += smarsa_BirthdaySpacings_nu(gen, res16, 100, 4 * MILLION / 10, 20, 16, 13, 1);
#endif
        sres_DeletePoisson(res16);

        return pVal;
    }
    case 17:
    {
        sres_Poisson *res17 = sres_CreatePoisson();
        std::string pVal = "";
#ifdef USE_LONGLONG
        pVal += smarsa_BirthdaySpacings_nu(gen, res17, 3, 20 * MILLION, 22, 256, 8, 1);
#else
        pVal += smarsa_BirthdaySpacings_nu(gen, res17, 100, 4 * MILLION / 10, 26, 16, 13, 1);
#endif
        sres_DeletePoisson(res17);

        return pVal;
    }
    case 18:
    {
        lebool flag = snpair_mNP2S_Flag;
        snpair_mNP2S_Flag = TRUE;
        snpair_Res *res18 = snpair_CreateRes();
        std::string pVal = snpair_ClosePairs_GetPVal_CPairs_nu(gen, res18, 10, 2 * MILLION, 0, 2, 0, 30, 10);
        snpair_DeleteRes(res18);
        snpair_mNP2S_Flag = flag;

        return pVal;
    }
    case 19:
    {
        lebool flag = snpair_mNP2S_Flag;
        snpair_mNP2S_Flag = TRUE;
        snpair_Res *res19 = snpair_CreateRes();
        std::string pVal = snpair_ClosePairs_GetPVal_CPairs_nu(gen, res19, 10, 2 * MILLION, 0, 3, 0, 30, 10);
        snpair_DeleteRes(res19);
        snpair_mNP2S_Flag = flag;

        return pVal;
    }
    case 20:
    {
        lebool flag = snpair_mNP2S_Flag;
        snpair_mNP2S_Flag = TRUE;
        snpair_Res *res20 = snpair_CreateRes();
        std::string pVal = snpair_ClosePairs_GetPVal_CPairs_nu(gen, res20, 5, 2 * MILLION, 0, 7, 0, 30, 10);
        snpair_DeleteRes(res20);
        snpair_mNP2S_Flag = flag;

        return pVal;
    }
    case 21:
    {
        lebool flag = snpair_mNP2S_Flag;
        snpair_mNP2S_Flag = TRUE;
        snpair_Res *res21 = snpair_CreateRes();
        std::string pVal = snpair_ClosePairsBitMatch_nu(gen, res21, 4, 4 * MILLION, 0, 2);
        snpair_DeleteRes(res21);
        snpair_mNP2S_Flag = flag;

        return pVal;
    }
    case 22:
    {
        lebool flag = snpair_mNP2S_Flag;
        snpair_mNP2S_Flag = TRUE;
        snpair_Res *res22 = snpair_CreateRes();
        std::string pVal = snpair_ClosePairsBitMatch_nu(gen, res22, 2, 4 * MILLION, 0, 4);
        snpair_DeleteRes(res22);
        snpair_mNP2S_Flag = flag;

        return pVal;
    }
    case 23:
    {
        sres_Chi2 *res23 = sres_CreateChi2();
        std::string pVal = sknuth_SimpPoker_nu(gen, res23, 1, 40 * MILLION, 0, 16, 16);
        sres_DeleteChi2(res23);

        return pVal;
    }
    case 24:
    {
        sres_Chi2 *res24 = sres_CreateChi2();
        std::string pVal = sknuth_SimpPoker_nu(gen, res24, 1, 40 * MILLION, 26, 16, 16);
        sres_DeleteChi2(res24);

        return pVal;
    }
    case 25:
    {
        sres_Chi2 *res25 = sres_CreateChi2();
        std::string pVal = sknuth_SimpPoker_nu(gen, res25, 1, 10 * MILLION, 0, 64, 64);
        sres_DeleteChi2(res25);

        return pVal;
    }
    case 26:
    {
        sres_Chi2 *res26 = sres_CreateChi2();
        std::string pVal = sknuth_SimpPoker_nu(gen, res26, 1, 10 * MILLION, 24, 64, 64);
        sres_DeleteChi2(res26);

        return pVal;
    }
    case 27:
    {
        sres_Chi2 *res27 = sres_CreateChi2();
        std::string pVal = sknuth_CouponCollector_nu(gen, res27, 1, 40 * MILLION, 0, 4);
        sres_DeleteChi2(res27);

        return pVal;
    }
    case 28:
    {
        sres_Chi2 *res28 = sres_CreateChi2();
        std::string pVal = sknuth_CouponCollector_nu(gen, res28, 1, 40 * MILLION, 28, 4);
        sres_DeleteChi2(res28);

        return pVal;
    }
    case 29:
    {
        sres_Chi2 *res29 = sres_CreateChi2();
        std::string pVal = sknuth_CouponCollector_nu(gen, res29, 1, 10 * MILLION, 0, 16);
        sres_DeleteChi2(res29);

        return pVal;
    }
    case 30:
    {
        sres_Chi2 *res30 = sres_CreateChi2();
        std::string pVal = sknuth_CouponCollector_nu(gen, res30, 1, 10 * MILLION, 26, 16);
        sres_DeleteChi2(res30);

        return pVal;
    }
    case 31:
    {
        sres_Chi2 *res31 = sres_CreateChi2();
        std::string pVal = sknuth_Gap_nu(gen, res31, 1, 100 * MILLION, 0, 0.0, 0.125);
        sres_DeleteChi2(res31);

        return pVal;
    }
    case 32:
    {
        sres_Chi2 *res32 = sres_CreateChi2();
        std::string pVal = sknuth_Gap_nu(gen, res32, 1, 100 * MILLION, 27, 0.0, 0.125);
        sres_DeleteChi2(res32);

        return pVal;
    }
    case 33:
    {
        sres_Chi2 *res33 = sres_CreateChi2();
        std::string pVal = sknuth_Gap_nu(gen, res33, 1, 5 * MILLION, 0, 0.0, 1.0 / 256.0);
        sres_DeleteChi2(res33);

        return pVal;
    }
    case 34:
    {
        sres_Chi2 *res34 = sres_CreateChi2();
        std::string pVal = sknuth_Gap_nu(gen, res34, 1, 5 * MILLION, 22, 0.0, 1.0 / 256.0);
        sres_DeleteChi2(res34);

        return pVal;
    }
    case 35:
    {
        sres_Chi2 *res35 = sres_CreateChi2();
        std::string pVal = sknuth_Run_nu_Mean(gen, res35, 1, 500 * MILLION, 0, TRUE);
        sres_DeleteChi2(res35);

        return pVal;
    }
    case 36:
    {
        sres_Chi2 *res36 = sres_CreateChi2();
        std::string pVal = sknuth_Run_nu_Mean(gen, res36, 1, 500 * MILLION, 15, FALSE);
        sres_DeleteChi2(res36);

        return pVal;
    }
    case 37:
    {
        sres_Chi2 *res37 = sres_CreateChi2();
        std::string pVal = sknuth_Permutation_nu(gen, res37, 1, 50 * MILLION, 0, 10);
        sres_DeleteChi2(res37);

        return pVal;
    }
    case 38:
    {
        sres_Chi2 *res38 = sres_CreateChi2();
        std::string pVal = sknuth_Permutation_nu(gen, res38, 1, 50 * MILLION, 15, 10);
        sres_DeleteChi2(res38);

        return pVal;
    }
    case 39:
    {
        sknuth_Res2 *res39 = sknuth_CreateRes2();
        std::string pVal = sknuth_CollisionPermut_nu(gen, res39, 5, 10 * MILLION, 0, 13);
        sknuth_DeleteRes2(res39);

        return pVal;
    }
    case 40:
    {
        sknuth_Res2 *res40 = sknuth_CreateRes2();
        std::string pVal = sknuth_CollisionPermut_nu(gen, res40, 5, 10 * MILLION, 15, 13);
        sknuth_DeleteRes2(res40);

        return pVal;
    }
    case 41:
    {
        sknuth_Res1 *res41 = sknuth_CreateRes1();
        std::string pVal = sknuth_MaxOft_nu_Sum_AD(gen, res41, 10, 10 * MILLION, 0, MILLION / 10, 5);
        sknuth_DeleteRes1(res41);

        return pVal;
    }
    case 42:
    {
        sknuth_Res1 *res42 = sknuth_CreateRes1();
        std::string pVal = sknuth_MaxOft_nu_Sum_AD(gen, res42, 5, 10 * MILLION, 0, MILLION / 10, 10);
        sknuth_DeleteRes1(res42);

        return pVal;
    }
    case 43:
    {
        sknuth_Res1 *res43 = sknuth_CreateRes1();
        std::string pVal = sknuth_MaxOft_nu_Mean_Mean(gen, res43, 1, 10 * MILLION, 0, MILLION / 10, 20);
        sknuth_DeleteRes1(res43);

        return pVal;
    }
    case 44:
    {
        sknuth_Res1 *res44 = sknuth_CreateRes1();
        std::string pVal = sknuth_MaxOft_nu_Mean_Mean(gen, res44, 1, 10 * MILLION, 0, MILLION / 10, 30);
        sknuth_DeleteRes1(res44);

        return pVal;
    }
    case 45:
    {
        sres_Basic *res45 = sres_CreateBasic();
        std::string pVal = svaria_SampleProd_nu_Mean(gen, res45, 1, 10 * MILLION, 0, 10);
        sres_DeleteBasic(res45);

        return pVal;
    }
    case 46:
    {
        sres_Basic *res46 = sres_CreateBasic();
        std::string pVal = svaria_SampleProd_nu_Mean(gen, res46, 1, 10 * MILLION, 0, 30);
        sres_DeleteBasic(res46);

        return pVal;
    }
    case 47:
    {
        sres_Basic *res47 = sres_CreateBasic();
        std::string pVal = svaria_SampleMean_nu(gen, res47, 10 * MILLION, 20, 0);
        sres_DeleteBasic(res47);

        return pVal;
    }
    case 48:
    {
        sres_Basic *res48 = sres_CreateBasic();
        std::string pVal = svaria_SampleCorr_nu(gen, res48, 1, 500 * MILLION, 0, 1);
        sres_DeleteBasic(res48);

        return pVal;
    }
    case 49:
    {
        sres_Basic *res49 = sres_CreateBasic();
        std::string pVal = svaria_AppearanceSpacings_nu(gen, res49, 1, 10 * MILLION, 400 * MILLION, r, 30, 15);
        sres_DeleteBasic(res49);

        return pVal;
    }
    case 50:
    {
        sres_Basic *res50 = sres_CreateBasic();
        std::string pVal = svaria_AppearanceSpacings_nu(gen, res50, 1, 10 * MILLION, 100 * MILLION, 20, 10, 15);
        sres_DeleteBasic(res50);

        return pVal;
    }
    case 51:
    {
        sres_Chi2 *res51 = sres_CreateChi2();
        std::string pVal = svaria_WeightDistrib_nu(gen, res51, 1, 2 * MILLION, 0, 256, 0.0, 0.125);
        sres_DeleteChi2(res51);

        return pVal;
    }
    case 52:
    {
        sres_Chi2 *res52 = sres_CreateChi2();
        std::string pVal = svaria_WeightDistrib_nu(gen, res52, 1, 2 * MILLION, 8, 256, 0.0, 0.125);
        sres_DeleteChi2(res52);

        return pVal;
    }
    case 53:
    {
        sres_Chi2 *res53 = sres_CreateChi2();
        std::string pVal = svaria_WeightDistrib_nu(gen, res53, 1, 2 * MILLION, 16, 256, 0.0, 0.125);
        sres_DeleteChi2(res53);

        return pVal;
    }
    case 54:
    {
        sres_Chi2 *res54 = sres_CreateChi2();
        std::string pVal = svaria_WeightDistrib_nu(gen, res54, 1, 2 * MILLION, 24, 256, 0.0, 0.125);
        sres_DeleteChi2(res54);

        return pVal;
    }
    case 55:
    {
        sres_Chi2 *res55 = sres_CreateChi2();
        std::string pVal = svaria_SumCollector_nu(gen, res55, 1, 20 * MILLION, 0, 10.0);
        sres_DeleteChi2(res55);

        return pVal;
    }
    case 56:
    {
        sres_Chi2 *res56 = sres_CreateChi2();
        std::string pVal = smarsa_MatrixRank_nu_Mean(gen, res56, 1, MILLION, r, s, 2 * s, 2 * s);
        sres_DeleteChi2(res56);

        return pVal;
    }
    case 57:
    {
        sres_Chi2 *res57 = sres_CreateChi2();
        std::string pVal = smarsa_MatrixRank_nu_Mean(gen, res57, 1, MILLION, 20, 10, 2 * s, 2 * s);
        sres_DeleteChi2(res57);

        return pVal;
    }
    case 58:
    {
        sres_Chi2 *res58 = sres_CreateChi2();
        std::string pVal = smarsa_MatrixRank_nu_Mean(gen, res58, 1, 50 * THOUSAND, r, s, 10 * s, 10 * s);
        sres_DeleteChi2(res58);

        return pVal;
    }
    case 59:
    {
        sres_Chi2 *res59 = sres_CreateChi2();
        std::string pVal = smarsa_MatrixRank_nu_Mean(gen, res59, 1, 50 * THOUSAND, 20, 10, 10 * s, 10 * s);
        sres_DeleteChi2(res59);

        return pVal;
    }
    case 60:
    {
        sres_Chi2 *res60 = sres_CreateChi2();
        std::string pVal = smarsa_MatrixRank_nu_Mean(gen, res60, 1, 2 * THOUSAND, r, s, 40 * s, 40 * s);
        sres_DeleteChi2(res60);

        return pVal;
    }
    case 61:
    {
        sres_Chi2 *res61 = sres_CreateChi2();
        std::string pVal = smarsa_MatrixRank_nu_Mean(gen, res61, 1, 2 * THOUSAND, 20, 10, 40 * s, 40 * s);
        sres_DeleteChi2(res61);

        return pVal;
    }
    case 62:
    {
        sres_Chi2 *res62 = sres_CreateChi2();
        std::string pVal = smarsa_Savir2_nu_Mean(gen, res62, 1, 20 * MILLION, 0, 1024 * 1024, 30);
        sres_DeleteChi2(res62);

        return pVal;
    }
    case 63:
    {
        smarsa_Res2 *res63 = smarsa_CreateRes2();
        std::string pVal = smarsa_GCD_nu_Mean(gen, res63, 1, 100 * MILLION, 0, 30);
        smarsa_DeleteRes2(res63);

        return pVal;
    }
    case 64:
    {
        smarsa_Res2 *res64 = smarsa_CreateRes2();
        std::string pVal = smarsa_GCD_nu_Mean(gen, res64, 1, 40 * MILLION, 10, 20);
        smarsa_DeleteRes2(res64);

        return pVal;
    }
    case 65:
    {
        swalk_Res *res65 = swalk_CreateRes();
        std::string pVal = swalk_RandomWalk1_nu(gen, res65, 1, 50 * MILLION, r, s, 90, 90, 1);
        swalk_DeleteRes(res65);

        return pVal;
    }
    case 66:
    {
        swalk_Res *res66 = swalk_CreateRes();
        std::string pVal = swalk_RandomWalk1_nu(gen, res66, 1, 10 * MILLION, 20, 10, 90, 90, 1);
        swalk_DeleteRes(res66);

        return pVal;
    }
    case 67:
    {
        swalk_Res *res67 = swalk_CreateRes();
        std::string pVal = swalk_RandomWalk1_nu(gen, res67, 1, 5 * MILLION, r, s, 1000, 1000, 1);
        swalk_DeleteRes(res67);

        return pVal;
    }
    case 68:
    {
        swalk_Res *res68 = swalk_CreateRes();
        std::string pVal = swalk_RandomWalk1_nu(gen, res68, 1, MILLION, 20, 10, 1000, 1000, 1);
        swalk_DeleteRes(res68);

        return pVal;
    }
    case 69:
    {
        swalk_Res *res69 = swalk_CreateRes();
        std::string pVal = swalk_RandomWalk1_nu(gen, res69, 1, MILLION / 2, r, s, 10000, 10000, 1);
        swalk_DeleteRes(res69);

        return pVal;
    }
    case 70:
    {
        swalk_Res *res70 = swalk_CreateRes();
        std::string pVal = swalk_RandomWalk1_nu(gen, res70, 1, MILLION / 10, 20, 10, 10000, 10000, 1);
        swalk_DeleteRes(res70);

        return pVal;
    }
    case 71:
    {
        scomp_Res *res71 = scomp_CreateRes();
        std::string pVal = scomp_LinearComp_nu(gen, res71, 1, 120 * THOUSAND, r, 1);
        scomp_DeleteRes(res71);

        return pVal;
    }
    case 72:
    {
        scomp_Res *res72 = scomp_CreateRes();
        std::string pVal = scomp_LinearComp_nu(gen, res72, 1, 120 * THOUSAND, 29, 1);
        scomp_DeleteRes(res72);

        return pVal;
    }
    case 73:
    {
        sres_Basic *res73 = sres_CreateBasic();
        std::string pVal = scomp_LempelZiv_nu(gen, res73, 10, 25, r, s);
        sres_DeleteBasic(res73);

        return pVal;
    }
    case 74:
    {
        sspectral_Res *res74 = sspectral_CreateRes();
        std::string pVal = sspectral_Fourier3_nu(gen, res74, 50 * THOUSAND, 14, r, s);
        sspectral_DeleteRes(res74);

        return pVal;
    }
    case 75:
    {
        sspectral_Res *res75 = sspectral_CreateRes();
        std::string pVal = sspectral_Fourier3_nu(gen, res75, 50 * THOUSAND, 14, 20, 10);
        sspectral_DeleteRes(res75);

        return pVal;
    }
    case 76:
    {
        sstring_Res2 *res76 = sstring_CreateRes2();
        std::string pVal = sstring_LongestHeadRun_nu(gen, res76, 1, 1000, r, s, 20 + 10 * MILLION);
        sstring_DeleteRes2(res76);

        return pVal;
    }
    case 77:
    {
        sstring_Res2 *res77 = sstring_CreateRes2();
        std::string pVal = sstring_LongestHeadRun_nu(gen, res77, 1, 300, 20, 10, 20 + 10 * MILLION);
        sstring_DeleteRes2(res77);

        return pVal;
    }
    case 78:
    {
        sres_Chi2 *res78 = sres_CreateChi2();
        std::string pVal = sstring_PeriodsInStrings_nu_Mean(gen, res78, 1, 300 * MILLION, r, s);
        sres_DeleteChi2(res78);

        return pVal;
    }
    case 79:
    {
        sres_Chi2 *res79 = sres_CreateChi2();
        std::string pVal = sstring_PeriodsInStrings_nu_Mean(gen, res79, 1, 300 * MILLION, 15, 15);
        sres_DeleteChi2(res79);

        return pVal;
    }
    case 80:
    {
        sres_Basic *res80 = sres_CreateBasic();
        std::string pVal = sstring_HammingWeight2_nu(gen, res80, 100, 100 * MILLION, r, s, MILLION);
        sres_DeleteBasic(res80);

        return pVal;
    }
    case 81:
    {
        sres_Basic *res81 = sres_CreateBasic();
        std::string pVal = sstring_HammingWeight2_nu(gen, res81, 30, 100 * MILLION, 20, 10, MILLION);
        sres_DeleteBasic(res81);

        return pVal;
    }
    case 82:
    {
        sstring_Res *res82 = sstring_CreateRes();
        std::string pVal = sstring_HammingCorr_nu(gen, res82, 1, 500 * MILLION, r, s, s);
        sstring_DeleteRes(res82);

        return pVal;
    }
    case 83:
    {
        sstring_Res *res83 = sstring_CreateRes();
        std::string pVal = sstring_HammingCorr_nu(gen, res83, 1, 50 * MILLION, r, s, 10 * s);
        sstring_DeleteRes(res83);

        return pVal;
    }
    case 84:
    {
        sstring_Res *res84 = sstring_CreateRes();
        std::string pVal = sstring_HammingCorr_nu(gen, res84, 1, 10 * MILLION, r, s, 40 * s);
        sstring_DeleteRes(res84);

        return pVal;
    }
    case 85:
    {
        sstring_Res *res85 = sstring_CreateRes();
        std::string pVal = sstring_HammingIndep_nu_Mean(gen, res85, 1, 300 * MILLION, r, s, s, 0);
        sstring_DeleteRes(res85);

        return pVal;
    }
    case 86:
    {
        sstring_Res *res86 = sstring_CreateRes();
        std::string pVal = sstring_HammingIndep_nu_Mean(gen, res86, 1, 100 * MILLION, 20, 10, s, 0);
        sstring_DeleteRes(res86);

        return pVal;
    }
    case 87:
    {
        sstring_Res *res87 = sstring_CreateRes();
        std::string pVal = sstring_HammingIndep_nu_Mean(gen, res87, 1, 30 * MILLION, r, s, 10 * s, 0);
        sstring_DeleteRes(res87);

        return pVal;
    }
    case 88:
    {
        sstring_Res *res88 = sstring_CreateRes();
        std::string pVal = sstring_HammingIndep_nu_Mean(gen, res88, 1, 10 * MILLION, 20, 10, 10 * s, 0);
        sstring_DeleteRes(res88);

        return pVal;
    }
    case 89:
    {
        sstring_Res *res89 = sstring_CreateRes();
        std::string pVal = sstring_HammingIndep_nu_Mean(gen, res89, 1, 10 * MILLION, r, s, 40 * s, 0);
        sstring_DeleteRes(res89);

        return pVal;
    }
    case 90:
    {
        sstring_Res *res90 = sstring_CreateRes();
        std::string pVal = sstring_HammingIndep_nu_Mean(gen, res90, 1, MILLION, 20, 10, 40 * s, 0);
        sstring_DeleteRes(res90);

        return pVal;
    }
    case 91:
    {
        sstring_Res3 *res91 = sstring_CreateRes3();
        std::string pVal = sstring_Run_nu(gen, res91, 1, 1 * BILLION, r, s);
        sstring_DeleteRes3(res91);

        return pVal;
    }
    case 92:
    {
        sstring_Res3 *res92 = sstring_CreateRes3();
        std::string pVal = sstring_Run_nu(gen, res92, 1, 1 * BILLION, 20, 10);
        sstring_DeleteRes3(res92);

        return pVal;
    }
    case 93:
    {
        sres_Basic *res93 = sres_CreateBasic();
        std::string pVal = sstring_AutoCor_nu(gen, res93, 10, 30 + BILLION, r, s, 1);
        sres_DeleteBasic(res93);

        return pVal;
    }
    case 94:
    {
        sres_Basic *res94 = sres_CreateBasic();
        std::string pVal = sstring_AutoCor_nu(gen, res94, 5, 1 + BILLION, 20, 10, 1);
        sres_DeleteBasic(res94);

        return pVal;
    }
    case 95:
    {
        sres_Basic *res95 = sres_CreateBasic();
        std::string pVal = sstring_AutoCor_nu(gen, res95, 10, 31 + BILLION, r, s, s);
        sres_DeleteBasic(res95);

        return pVal;
    }
    case 96:
    {
        sres_Basic *res96 = sres_CreateBasic();
        std::string pVal = sstring_AutoCor_nu(gen, res96, 5, 11 + BILLION, 20, 10, 10);
        sres_DeleteBasic(res96);

        return pVal;
    }
    default:
    {
        std::string pVal = "Invalid case.\n";

        return pVal;
    }
    }
}
