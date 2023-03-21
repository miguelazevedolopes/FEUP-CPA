#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <string>
#include <bits/stdc++.h>
#include <omp.h>
#include <papi.h>

#define N_THREADS 8

int ret;
int EventSet = PAPI_NULL;

void SieveOfEratosthenes(long long n, std::ofstream &outputFile)
{
    // Array in style -> [2, 3, 4, 5, 6, 7, 8, ..., n]
    // Size is n - 1
    long long size = n - 1;
    bool *marks = new bool[size];

    std::cout << "Calculating SieveOfEratosthenes" << std::endl;

    // Start time
    auto t1 = std::chrono::high_resolution_clock::now();

    ret = PAPI_start(EventSet);
    if (ret != PAPI_OK)
      std::cout << "ERROR: Start PAPI" << std::endl;
    // Calculate primes
    for (long long k = 0; (k + 2) * (k + 2) < n; k++)
    {
        // Number is always 2 + index of array
        long long num = k + 2;

        // Go from num * num until n to check for multiples
        for (long long j = num * num - 2; j < size; j++)
        {
            // Number is always 2 + index of array
            long long multNum = j + 2;

            // Mark number if divisible by study number
            marks[j] |= (multNum % num == 0);
        }
    }
    long long values[3];
    ret = PAPI_stop(EventSet, values);
    if (ret != PAPI_OK)
        std::cout << "ERROR: Stop PAPI" << std::endl;

    ret = PAPI_reset(EventSet);
    if (ret != PAPI_OK)
        std::cout << "FAIL reset" << std::endl;

    auto t2 = std::chrono::high_resolution_clock::now();

    // Output Execution time
    outputFile << "Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0 << "s" << std::endl << "Total Instructions: " << values[0] << std::endl << "L1 Cache Misses: " << values[1] << std::endl << "L2 Cache Misses: " << values[2] << std::endl
               << std::endl;
    std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0 << "s" << std::endl
              << std::endl;

    // Output Primes
    for (long long k = 0; k < size; k++)
    {
        // Unmarked are prime
        if (!marks[k])
        {
            // Number is always 2 + index of array
            outputFile << k + 2 << std::endl;
        }
    }

    // Clear memory :D
    delete marks;
}

void SieveOfEratosthenesFastMarking(long long n, std::ofstream &outputFile)
{
    // Array in style -> [2, 3, 4, 5, 6, 7, 8, ..., n]
    // Size is n - 1
    long long size = n - 1;
    bool *marks = new bool[size];

    std::cout << "Calculating SieveOfEratosthenesFastMarking" << std::endl;

    // Start time
    auto t1 = std::chrono::high_resolution_clock::now();

    ret = PAPI_start(EventSet);
    if (ret != PAPI_OK)
      std::cout << "ERROR: Start PAPI" << std::endl;
    // Calculate primes
    for (long long k = 0; (k + 2) * (k + 2) < n; k++)
    {
        // Number is always 2 + index of array
        long long num = k + 2;

        // Go from num * num until n to check for multiples
        long long counter = 1;
        for (long long multNum = num * num; multNum <= n; counter++)
        {
            // Mark number if divisible by study number
            marks[multNum - 2] |= (multNum % num == 0);

            // Next number is always the base num * num plus the multiple of the base number
            multNum = num * num + num * counter;
        }
    }

    long long values[3];
    ret = PAPI_stop(EventSet, values);
    if (ret != PAPI_OK)
        std::cout << "ERROR: Stop PAPI" << std::endl;

    ret = PAPI_reset(EventSet);
    if (ret != PAPI_OK)
        std::cout << "FAIL reset" << std::endl;

    auto t2 = std::chrono::high_resolution_clock::now();

    // Output Execution time
    outputFile << "Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0 << "s" << std::endl << "Total Instructions: " << values[0] << std::endl << "L1 Cache Misses: " << values[1] << std::endl << "L2 Cache Misses: " << values[2] << std::endl
               << std::endl;
    std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0 << "s" << std::endl
              << std::endl;

    // Output Primes
    for (long long k = 0; k < size; k++)
    {
        // Unmarked are prime
        if (!marks[k])
        {
            // Number is always 2 + index of array
            outputFile << k + 2 << std::endl;
        }
    }



    // Clear memory :D
    delete marks;
}

void SieveOfEratosthenesBlock(long long n, std::ofstream &outputFile)
{
    std::cout << "Calculating SieveOfEratosthenesBlock" << std::endl;

    // Array in style -> [2, 3, 4, 5, 6, 7, 8, ..., n]
    // Size is n - 1
    long long limit = floor(sqrt(n)) + 1;

    // First set of primes, used to check for primes in other blocks
    long long* prime = new long long[limit];
    long long primeSize = 0;

    bool* primeFull = new bool[n];
    memset (primeFull, true, sizeof (bool) * n);


    bool marks[limit+1];
    memset(marks, false, sizeof(marks));

    // Start time
    auto t1 = std::chrono::high_resolution_clock::now();

    ret = PAPI_start(EventSet);
    if (ret != PAPI_OK)
      std::cout << "ERROR: Start PAPI" << std::endl;
    for (long long p=2; p*p<limit; p++)
    {
        // If p is not changed, then it is a prime
        if (marks[p] == false)
        {
            // Update all multiples of p
            for (long long i=p*p; i<limit; i+=p)
                marks[i] = true;
        }
    }
 

    // Print all prime numbers and store them in prime
    for (long long p=2; p<limit; p++)
    {
        if (!marks[p])
        {
            prime[primeSize]=p;
            primeSize++;
        }
    }


    long long low = limit;
    long long high = 2*limit;

    for (low,high; low < n;low+=limit,high+=limit)
    {
        if (high >= n)
           high = n;

        bool mark[limit+1];
        memset(mark, false, sizeof(mark));


        for (long long i = 0; i < primeSize; i++)
        {

            long long loLim = floor(low/prime[i]) * prime[i];
            if (loLim < low)
                loLim += prime[i];
 

            for (long long j=loLim; j<high; j+=prime[i])
                mark[j-low] = true;
        }

        for (long long i = low; i<high; i++){
            if (mark[i - low] == false){
                primeFull[i]=false;
            }
        }
                
    }
    long long values[3];
    ret = PAPI_stop(EventSet, values);
    if (ret != PAPI_OK)
        std::cout << "ERROR: Stop PAPI" << std::endl;

    ret = PAPI_reset(EventSet);
    if (ret != PAPI_OK)
        std::cout << "FAIL reset" << std::endl;

    auto t2 = std::chrono::high_resolution_clock::now();

    // Output Execution time
    outputFile << "Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0 << "s" << std::endl << "Total Instructions: " << values[0] << std::endl << "L1 Cache Misses: " << values[1] << std::endl << "L2 Cache Misses: " << values[2] << std::endl
               << std::endl;
    std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0 << "s" << std::endl
              << std::endl;

    for (long long i=0;i<primeSize;i++){
        outputFile << prime[i]<< std::endl;
    }
    for (long long i=0;i<n;i++){
        if(!primeFull[i]){
            outputFile << i << std::endl;
        }
    }
    delete prime,primeFull;


}


// Applying the SPMD model using OpenMP 
void SieveOfEratosthenesBlockOMP(long long n, std::ofstream &outputFile)
{
    std::cout << "Calculating SieveOfEratosthenesBlockOMP" << std::endl;

    // Array in style -> [2, 3, 4, 5, 6, 7, 8, ..., n]
    // Size is n - 1
    long long limit = floor(sqrt(n)) + 1;

    // First set of primes, used to check for primes in other blocks
    long long* prime = new long long[limit];
    long long primeSize = 0;

    bool* primeFull = new bool[n];
    memset (primeFull, true, sizeof (bool) * n);

    bool marks[limit+1];
    memset(marks, false, sizeof(marks));

    // Start time
    auto t1 = std::chrono::high_resolution_clock::now();

    for (long long p=2; p*p<limit; p++)
    {
        // If p is not changed, then it is a prime
        if (marks[p] == false)
        {
            // Update all multiples of p
            for (long long i=p*p; i<limit; i+=p)
                marks[i] = true;
        }
    }

    // Print all prime numbers and store them in prime
    for (long long p=2; p<limit; p++)
    {
        if (!marks[p])
        {
            prime[primeSize]=p;
            primeSize++;
        }
    }


    omp_set_num_threads(N_THREADS);
    #pragma omp parallel shared(primeFull)
    {
        int id, thread_limit,thread_n;
        id = omp_get_thread_num();

        // 40 000 is the size of each segment handled by a thread/core.
        // Considering 1 bool occupies about 1 byte and L1 cache in my CPU is 64 kilobytes (64 000 bytes)
        // 40 000 was the threshold value that presented better results
        // (higher and lower values had worst results )
        thread_limit =  limit<40000 ? limit : 40000; 

        long long low = (long long)floor(n/N_THREADS)*(id)+limit;
        long long high = low + thread_limit;
        thread_n = id==(N_THREADS-1) ? n : (long long)floor(n/N_THREADS)*(id+1)+limit;

        for (low,high; low < thread_n;low+=thread_limit,high+=thread_limit)
        {
            if (high >= thread_n)
                high = thread_n;

            bool mark[thread_limit+1];
            memset(mark, false, sizeof(mark));

            for (long long i = 0; i < primeSize; i++)
            {
                long long loLim = floor(low/prime[i]) * prime[i];
                if (loLim < low)
                    loLim += prime[i];


                for (long long j=loLim; j<high; j+=prime[i])
                    mark[j-low] = true;

            }

            for (long long i = low; i<high; i++){
                if (mark[i - low] == false){
                    primeFull[i]=false;
                }
            }
        }                  
    }
    

    auto t2 = std::chrono::high_resolution_clock::now();

    // Output Execution time
    outputFile << "Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0 << "s" << std::endl
               << std::endl;
    std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0 << "s" << std::endl
              << std::endl;

    for (long long i=0;i<primeSize;i++){
        outputFile << prime[i]<< std::endl;
    }
    for (long long i=0;i<n;i++){

        if(!primeFull[i]){
            outputFile << i << std::endl;
        }
    }

    delete prime;
    delete primeFull;
}

void SieveOfEratosthenesBlockOMPTask(long long n, std::ofstream &outputFile)
{
    std::cout << "Calculating SieveOfEratosthenesBlockOMPTask" << std::endl;

    // Array in style -> [2, 3, 4, 5, 6, 7, 8, ..., n]
    // Size is n - 1
    long long limit = floor(sqrt(n)) + 1;

    // First set of primes, used to check for primes in other blocks
    long long* prime = new long long[limit];
    long long primeSize = 0;

    bool* primeFull = new bool[n];
    memset (primeFull, true, sizeof (bool) * n);


    bool marks[limit+1];
    memset(marks, false, sizeof(marks));

    // Start time
    auto t1 = std::chrono::high_resolution_clock::now();

    long long stop_limit =floor(sqrt(limit)) + 1;
    #pragma omp parallel
    #pragma omp single
    #pragma omp taskloop num_tasks(20)
        for (long long p=2; p<stop_limit; p++)
        {
            // If p is not changed, then it is a prime
            if (marks[p] == false)
            {
                // Update all multiples of p
                for (long long i=p*p; i<limit; i+=p)
                    marks[i] = true;
            }
        }
 

    // Print all prime numbers and store them in prime
    for (long long p=2; p<limit; p++)
    {
        if (!marks[p])
        {
            prime[primeSize]=p;
            primeSize++;
        }
    }



    #pragma omp parallel
    #pragma omp single
    #pragma omp taskloop num_tasks(20)
        for (long long low=limit; low < n;low+=limit){
            long long high=low+limit;
            if (high >= n)
                high = n;

            bool mark[limit+1];
            memset(mark, false, sizeof(bool)* (limit+1));

            for (long long i = 0; i < primeSize; i++)
            {

                long long loLim = floor(low/prime[i]) * prime[i];

                if (loLim < low)
                    loLim += prime[i];



                for (long long j=loLim; j<high; j+=prime[i]){

                    mark[j-low] = true;
                }
                    

            }

            for (long long i = low; i<high; i++){
                if (mark[i - low] == false){
                    primeFull[i]=false;
                }
            }
            
        }


    auto t2 = std::chrono::high_resolution_clock::now();

    // Output Execution time
    outputFile << "Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0 << "s" << std::endl
               << std::endl;
    std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0 << "s" << std::endl
              << std::endl;

    for (long long i=0;i<primeSize;i++){
        outputFile << prime[i]<< std::endl;
    }
    for (long long i=0;i<n;i++){
        if(!primeFull[i]){
            outputFile << i << std::endl;
        }
    }
    delete prime,primeFull;
}

void value_testing() {
    long long n;
    std::cout << "N: ";
    std::cin >> n;

    std::ofstream outputFile;
    std::cout << n << std::endl;

    outputFile.open("simple" + std::to_string(n) + ".txt");
    SieveOfEratosthenes(n, outputFile);
    outputFile.close();

    outputFile.open("fastMarking" + std::to_string(n) + ".txt");
    SieveOfEratosthenesFastMarking(n, outputFile);
    outputFile.close();

    outputFile.open("block" + std::to_string(n) + ".txt");
    SieveOfEratosthenesBlock(n, outputFile);
    outputFile.close();

    outputFile.open("blockOpm"+ std::to_string(n) + ".txt");
    SieveOfEratosthenesBlockOMP(n, outputFile);
    outputFile.close();

    outputFile.open("blockOpmTask"+ std::to_string(n) + ".txt");
    SieveOfEratosthenesBlockOMPTask(n, outputFile);
    outputFile.close();
}

void main_testing() {
    long long n;
    std::ofstream outputFile;
    for (int i = 25; i <= 32; i++) {
        n = pow(2, i);
        std::cout << n << std::endl;

        outputFile.open("simple10pow" + std::to_string(i) + ".txt");
        SieveOfEratosthenes(n, outputFile);
        outputFile.close();

        outputFile.open("fastMarking10pow" + std::to_string(i) + ".txt");
        SieveOfEratosthenesFastMarking(n, outputFile);
        outputFile.close();

        outputFile.open("block10pow" + std::to_string(i) + ".txt");
        SieveOfEratosthenesBlock(n, outputFile);
        outputFile.close();
    }
}

int main(int argc, char *argv[])
{
    ret = PAPI_library_init(PAPI_VER_CURRENT);
    if (ret != PAPI_VER_CURRENT)
        std::cout << "FAIL init" << std::endl;

    ret = PAPI_create_eventset(&EventSet);
    if (ret != PAPI_OK)
        std::cout << "ERROR: create eventset" << std::endl;

    ret = PAPI_add_event(EventSet, PAPI_TOT_INS);
    if (ret != PAPI_OK)
        std::cout << "ERROR: PAPI_TOT_INS" << std::endl;

    ret = PAPI_add_event(EventSet, PAPI_L1_DCM);
    if (ret != PAPI_OK)
        std::cout << "ERROR: PAPI_L1_DCM" << std::endl;

    ret = PAPI_add_event(EventSet, PAPI_L2_DCM);
    if (ret != PAPI_OK)
        std::cout << "ERROR: PAPI_L2_DCM" << std::endl;

    int i;
    std::cout << "Mode\n(1) - Value Testing\n(2) - Main Testing\n:";
    std::cin >> i;
    switch (i)
    {
    case 1:
        value_testing();
        break;
    case 2:
        main_testing();
    default:
        std::cout << "Invalid Value Aborting..." << std::endl;
        break;
    }

    ret = PAPI_remove_event(EventSet, PAPI_TOT_INS);
    if (ret != PAPI_OK)
        std::cout << "FAIL remove event PAPI_TOT_INS" << std::endl;

    ret = PAPI_remove_event(EventSet, PAPI_L1_DCM);
    if (ret != PAPI_OK)
        std::cout << "FAIL remove event PAPI_L1_DCM" << std::endl;

    ret = PAPI_remove_event(EventSet, PAPI_L2_DCM);
    if (ret != PAPI_OK)
        std::cout << "FAIL remove event PAPI_L2_DCM" << std::endl;

    ret = PAPI_destroy_eventset(&EventSet);
    if (ret != PAPI_OK)
        std::cout << "FAIL destroy" << std::endl;
}