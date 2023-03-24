#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <string>
#include <bits/stdc++.h>
#include <omp.h>
#include <papi.h>


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

    // Divide the range [0..n-1] in segments sized as sqrt(n)
    long long limit = floor(sqrt(n)) + 1;

    // First set of primes, used to check for primes in other blocks
    long long* prime = new long long[limit];
    long long primeSize = 0;

    // Array that holds primes 
    bool* primeFull = new bool[n];
    memset (primeFull, false, sizeof (bool) * n);


    // Start time
    auto t1 = std::chrono::high_resolution_clock::now();

    ret = PAPI_start(EventSet);
    if (ret != PAPI_OK)
      std::cout << "ERROR: Start PAPI" << std::endl;

    // Starts calculating the first set of primes
    for (long long p=2; p*p<limit; p++)
    {
        // If p is not changed, then it is a prime
        if (primeFull[p] == false)
        {
            // Update all multiples of p
            for (long long i=p*p; i<limit; i+=p)
                primeFull[i] = true;
        }
    }

    // Fill the prime array with the prime values (instead of just booleans in an array)
    for (long long p=2; p<limit; p++)
    {
        if (!primeFull[p])
        {
            prime[primeSize]=p;
            primeSize++;
        }
    }


    long long low = limit;
    long long high = 2*limit;

    // Start calculating the primes in range [limit,n-1]
    for (low,high; low < n;low+=limit,high+=limit)
    {
        if (high >= n)
           high = n;


        for (long long i = 0; i < primeSize; i++)
        {
            // Find the minimum number in range [low,high] that is
            // a multiple of prime[i]
            // For example, if low is 11 and prime[i] is 3,
            // we start with 12, so loLim=12.
            long long loLim = floor(low/prime[i]) * prime[i];
            if (loLim < low)
                loLim += prime[i];
 
            // Next starting with the loLim we mark every number
            // in range [loLim,high] that is a multiple of prime[i]
            // Using the previous example, it would mark 12, 15, 18,...
            // until it reaches high
            for (long long j=loLim; j<high; j+=prime[i])
                primeFull[j] = true;
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

    for (long long i=2;i<n;i++){

        if(!primeFull[i]){
            outputFile << i << std::endl;
        }
    }
    delete prime,primeFull;


}


// Applying the SPMD model using OpenMP 
void SieveOfEratosthenesBlockOMP(long long n,int n_cores, std::ofstream &outputFile)
{
    std::cout << "Calculating SieveOfEratosthenesBlockOMP" << std::endl;

    // Divide the range [0..n-1] in segments sized as sqrt(n)
    long long limit = floor(sqrt(n)) + 1;

    // First set of primes, used to check for primes in other blocks
    long long* prime = new long long[limit];
    long long primeSize = 0;

    // Array that holds primes 
    bool* primeFull = new bool[n];
    memset (primeFull, false, sizeof (bool) * n);


    // Start time
    auto t1 = std::chrono::high_resolution_clock::now();

    ret = PAPI_start(EventSet);
    if (ret != PAPI_OK)
      std::cout << "ERROR: Start PAPI" << std::endl;

    // Starts calculating the first set of primes
    for (long long p=2; p*p<limit; p++)
    {
        // If p is not changed, then it is a prime
        if (primeFull[p] == false)
        {
            // Update all multiples of p
            for (long long i=p*p; i<limit; i+=p)
                primeFull[i] = true;
        }
    }

    // Fill the prime array with the prime values (instead of just booleans in an array)
    for (long long p=2; p<limit; p++)
    {
        if (!primeFull[p])
        {
            prime[primeSize]=p;
            primeSize++;
        }
    }

    // Start calculating the primes in range [limit,n-1]
    omp_set_num_threads(n_cores);
    #pragma omp parallel shared(primeFull)
    {
        long long id,thread_n;
        id = omp_get_thread_num();


        // In each iteration the program will calculate prime numbers
        // in range [low,high]. These values are updated after each
        // iteration
        long long low = (long long)floor(n/n_cores)*(id)+limit;
        long long high = low + limit;
        thread_n = id==(n_cores-1) ? n : (long long)floor(n/n_cores)*(id+1)+limit;

        for (low,high; low < thread_n;low+=limit,high+=limit)
        {
            if (high >= thread_n)
                high = thread_n;


            for (long long i = 0; i < primeSize; i++)
            {
                // Find the minimum number in range [low,high] that is
                // a multiple of prime[i]
                // For example, if low is 11 and prime[i] is 3,
                // we start with 12, so loLim=12.
                long long loLim = floor(low/prime[i]) * prime[i];
                if (loLim < low)
                    loLim += prime[i];

                // Next starting with the loLim we mark every number
                // in range [loLim,high] that is a multiple of prime[i]
                // Using the previous example, it would mark 12, 15, 18,...
                // until it reaches high
                for (long long j=loLim; j<high; j+=prime[i]){
                    primeFull[j]=true;
                }

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

    
    for (long long i=2;i<n;i++){

        if(!primeFull[i]){
            outputFile << i << std::endl;
        }
    }

    delete prime;
    delete primeFull;
}

void SieveOfEratosthenesBlockOMPTask(long long n,int n_tasks, std::ofstream &outputFile)
{
    std::cout << "Calculating SieveOfEratosthenesBlockOMPTask" << std::endl;

    // Divide the range [0..n-1] in segments sized as sqrt(n)
    long long limit = floor(sqrt(n)) + 1;

    // First set of primes, used to check for primes in other blocks
    long long* prime = new long long[limit];
    long long primeSize = 0;

    // Array that holds primes 
    bool* primeFull = new bool[n];
    memset (primeFull, false, sizeof (bool) * n);


    // Start time
    auto t1 = std::chrono::high_resolution_clock::now();

    ret = PAPI_start(EventSet);
    if (ret != PAPI_OK)
      std::cout << "ERROR: Start PAPI" << std::endl;

    // Starts calculating the first set of primes
    for (long long p=2; p*p<limit; p++)
    {
        // If p is not changed, then it is a prime
        if (primeFull[p] == false)
        {
            // Update all multiples of p
            for (long long i=p*p; i<limit; i+=p)
                primeFull[i] = true;
        }
    }

    // Fill the prime array with the prime values (instead of just booleans in an array)
    for (long long p=2; p<limit; p++)
    {
        if (!primeFull[p])
        {
            prime[primeSize]=p;
            primeSize++;
        }
    }


    // Start calculating the primes in range [limit,n-1]
    #pragma omp parallel
    #pragma omp single
    #pragma omp taskloop num_tasks(n_tasks)
        for (long long low=limit; low < n;low+=limit){
            long long high=low+limit;
            if (high >= n)
                high = n;


            for (long long i = 0; i < primeSize; i++)
            {
                // Find the minimum number in range [low,high] that is
                // a multiple of prime[i]
                // For example, if low is 11 and prime[i] is 3,
                // we start with 12, so loLim=12.
                long long loLim = floor(low/prime[i]) * prime[i];
                if (loLim < low)
                    loLim += prime[i];

                // Next starting with the loLim we mark every number
                // in range [loLim,high] that is a multiple of prime[i]
                // Using the previous example, it would mark 12, 15, 18,...
                // until it reaches high
                for (long long j=loLim; j<high; j+=prime[i]){
                    primeFull[j] = true;
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

    for (long long i=2;i<n;i++){

        if(!primeFull[i]){
            outputFile << i << std::endl;
        }
    }

    delete prime;
    delete primeFull;
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
    SieveOfEratosthenesBlockOMP(n,8,outputFile);
    outputFile.close();

    outputFile.open("blockOpmTask"+ std::to_string(n) + ".txt");
    SieveOfEratosthenesBlockOMPTask(n,8,outputFile);
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
        
        if (i == 32) {
            for (int p = 1; p <= omp_get_num_procs(); p++) {
                outputFile.open("blockOpm10pow"+ std::to_string(i) + "P" + std::to_string(p) + ".txt");
                SieveOfEratosthenesBlockOMP(n, p, outputFile);
                outputFile.close();
            }

            for (int t = 1; t <= 50; t < 12 ? t++ : (t == 12 ? t += 8: t += 10)) {
                outputFile.open("blockOpmTask10pow"+ std::to_string(i) + "T" + std::to_string(t) + ".txt");
                SieveOfEratosthenesBlockOMPTask(n, t, outputFile);
                outputFile.close();
            }
        }
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
        break;
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