#include <mpi.h>
#include <iostream>
#include <vector>
#include <fstream>
#include "crushselect.h"

enum CRUSH_TYPE
{
    SMALL_CRUSH,
    CRUSH,
    BIG_CRUSH
};

void send_string(const std::string &str, int dest, int tag, MPI_Comm comm)
{
    int length = str.size();
    MPI_Send(&length, 1, MPI_INT, dest, tag, comm);
    MPI_Send(str.c_str(), length, MPI_CHAR, dest, tag + 1, comm);
}

std::string recv_string(int source, int tag, MPI_Comm comm)
{
    int length;
    MPI_Recv(&length, 1, MPI_INT, source, tag, comm, MPI_STATUS_IGNORE);

    std::string str(length, '\0');
    MPI_Recv(&str[0], length, MPI_CHAR, source, tag + 1, comm, MPI_STATUS_IGNORE);

    return str;
}

std::string MPI_Crush_Results(CRUSH_TYPE ct)
{
    int length = 0;
    switch (ct)
    {
    case SMALL_CRUSH:
        length = 11;
        break;
    case CRUSH:
        length = 97;
        break;
    case BIG_CRUSH:
        length = 107;
        break;
    default:
        std::cout << "Rank 0: Invalid test type." << std::endl;
        return std::string("");
    }
    auto start = std::chrono::high_resolution_clock::now();
    std::string pVals = "";
    for (int i = 1; i < length; ++i)
    {
        pVals += recv_string(i, 0, MPI_COMM_WORLD);
        pVals += ",";
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double>(end - start);
    pVals += std::to_string(duration.count()) + "\n";
    return pVals;
}

void MPI_Crush_Test(CRUSH_TYPE ct, unif01_Gen *gen, int testNum)
{
    std::string pVal = "";
    switch (ct)
    {
    case SMALL_CRUSH:
        pVal = selectSmallCrushTest(gen, testNum);
        send_string(pVal, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&testNum, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
        break;
    case CRUSH:
        pVal = selectCrushTest(gen, testNum);
        send_string(pVal, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&testNum, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
        break;
    case BIG_CRUSH:
        pVal = selectBigCrushTest(gen, testNum);
        send_string(pVal, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&testNum, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
        break;
    default:
        std::cout << "Test " << testNum << ": Invalid test type." << std::endl;
        break;
    }
}

void run_master(int world_size, CRUSH_TYPE ct)
{
    std::vector<std::string> results;
    bool logs = false;
    std::ofstream MasterFile("master.out", std::ios::out);

    int work = 0;
    switch (ct)
    {
    case SMALL_CRUSH:
        work = 10;
        break;
    case CRUSH:
        work = 96;
        break;
    case BIG_CRUSH:
        work = 106;
        break;
    default:
        if(logs)
            MasterFile << "[Master] Invalid test type." << std::endl;
        return;
    }

    int completed_work = 0;
    results.resize(work);

    if(logs)
        MasterFile << "[Master] Starting and distributing tasks..." << std::endl;

    // Initial task distribution
    int i;
    for (i = 1; i < world_size && i < work + 1; ++i)
    {
        MPI_Send(&i, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        if(logs)
            MasterFile << "[Master] Worker " << i << " sent test " << i << std::endl;
    }
    int work_assigned = i - 1;

    auto start = std::chrono::high_resolution_clock::now();
    std::string pVals = "";

    if(logs)
        MasterFile << "[Master] Beginning polling loop..." << std::endl;
    while (completed_work < work)
    {
        int flag = 0;
        MPI_Status status;

        // Check for any incoming messages
        MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);

        if (flag)
        {
            std::string result = recv_string(status.MPI_SOURCE, 0, MPI_COMM_WORLD);
            int testNum;
            MPI_Recv(&testNum, 1, MPI_INT, status.MPI_SOURCE, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            results[testNum - 1] = result;

            if(logs)
                MasterFile << "[Master] Worker " << status.MPI_SOURCE << " reported result: " << result << std::endl;
            completed_work++;

            if (work_assigned < work)
            {
                ++work_assigned;
                MPI_Send(&work_assigned, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
                if(logs)
                    MasterFile << "[Master] Worker " << status.MPI_SOURCE << " sent test " << work_assigned << std::endl;
            }
            else
                usleep(1);
            if(logs)
                MasterFile << "[Master] Work Complete: " << completed_work << " | Work Assigned: " << work_assigned << std::endl;
        }
    }
    if(logs)
        MasterFile << "[Master] All work done. Sending shutdown signals." << std::endl;

    for (int i = 1; i < world_size; ++i)
    {
        int die_signal = -1; // Using -1 as a special "Poison Pill" tag or value
        MPI_Send(&die_signal, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double>(end - start);
    for (auto r : results)
    {
        std::cout << r << ',';
    }
    std::cout << duration.count() << std::endl;
    
    MasterFile.close();
}

void run_worker(int rank, CRUSH_TYPE ct, unif01_Gen *gen)
{
    while (true)
    {
        int task;
        MPI_Recv(&task, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (task == -1)
        {
            break;
        }
        MPI_Crush_Test(ct, gen, task);
        
    }
}
