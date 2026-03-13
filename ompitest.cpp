#include "sprng_cpp.h"
#include "crushmpi.h"
#include <random>

Sprng *ptr;

/* Returns generated double-precision numbers for testU01 generator */
double genRanDbl(void)
{
    return ptr->sprng();
}

int numberofparameters(int gtype)
{
    int numparams = 0;
    if (gtype == 3 or gtype == 2)
    {
        numparams = 3;
    }
    else if (gtype == 1)
    {
        numparams = 7;
    }
    else if (gtype == 0 or gtype == 4)
    {
        numparams = 11;
    }
    else
    {
        numparams = 1;
    }
    return numparams;
}

string nameofGen(int gtype)
{
    string Name = "";
    if (gtype == 0)
    {
        Name = "LFG";
    }
    else if (gtype == 1)
    {
        Name = "LCG";
    }
    else if (gtype == 2)
    {
        Name = "LCG64";
    }
    else if (gtype == 3)
    {
        Name = "CMRG";
    }
    else if (gtype == 4)
    {
        Name = "MLFG";
    }
    else if (gtype == 5)
    {
        Name = "PMLCG";
    }
    else
    {
        Name = "NonSPRNGGen";
    }
    return Name;
}

int main(int argc, char *argv[])
{
    close(2);
    swrite_Basic = FALSE;
    swrite_Host = FALSE;
    swrite_Classes = FALSE;
    swrite_Parameters = FALSE;
    swrite_Collectors = FALSE;
    swrite_Counters = FALSE;

    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    int testType;

    if (argc < 3)
    {
        std::cout << "Usage: srun -n N " << argv[0] << " testType[0-2] genType[0-5]\n";
        return 1;
    }

    testType = atoi(argv[1]);

    if (testType < 0 || testType > 2)
    {
        std::cout << "Invalid test type." << std::endl;
        return 1;
    }

    if (world_rank == 0)
    {
        std::random_device rd;
        unsigned int seed = rd();
        std::mt19937 gen(seed);
        for (int i = 1; i < world_size; ++i)
        {
            int randNum = gen();
            MPI_Send(&randNum, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
        run_master(world_size, static_cast<CRUSH_TYPE>(testType));

        // std::cout << MPI_Crush_Results(static_cast<CRUSH_TYPE>(testType));
        // std::cout << "base seed: " << seed << std::endl;
    }
    else
    {

        unif01_Gen *gen;
        int streamnum, nstreams, gtype, seed, param;

        gtype = atoi(argv[2]);
        MPI_Recv(&seed, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        param = SPRNG_DEFAULT;
        streamnum = 0;
        nstreams = 1;

        ptr = SelectType(gtype);
        ptr->init_sprng(streamnum, nstreams, seed, param);

        char genName[] = "SPRNG Double";

        gen = unif01_CreateExternGen01(genName, genRanDbl);

        run_worker(world_rank, static_cast<CRUSH_TYPE>(testType), gen);

        // MPI_Crush_Test(static_cast<CRUSH_TYPE>(testType), gen, world_rank);

        // std::cout << "seed " << world_rank << ": " << seed << std::endl;
    }

    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
}
