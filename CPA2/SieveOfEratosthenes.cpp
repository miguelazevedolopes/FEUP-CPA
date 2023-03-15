#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <vector>
#include <bits/stdc++.h>


void SieveOfEratosthenes(long long n, std::ofstream &outputFile)
{
    // Array in style -> [2, 3, 4, 5, 6, 7, 8, ..., n]
    // Size is n - 1
    long long size = n - 1;
    bool *marks = new bool[size];

    std::cout << "Calculating SieveOfEratosthenes" << std::endl;

    // Start time
    auto t1 = std::chrono::high_resolution_clock::now();

    // Calculate primes
    for (long long k = 0; k < size; k++)
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

    auto t2 = std::chrono::high_resolution_clock::now();

    // Output Execution time
    outputFile << "Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0 << "s" << std::endl
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

    // Calculate primes
    for (long long k = 0; k < size; k++)
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

    auto t2 = std::chrono::high_resolution_clock::now();

    // Output Execution time
    outputFile << "Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0 << "s" << std::endl
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

// TODO
void SieveOfEratosthenesFastMarkingReorganized(long long n, std::ofstream &outputFile)
{

    std::cout << "Calculating SieveOfEratosthenesFastMarking" << std::endl;



    // Array in style -> [2, 3, 4, 5, 6, 7, 8, ..., n]
    // Size is n - 1
    long long limit = floor(sqrt(n)) + 1;

    std::vector<long long> prime;
    prime.reserve(limit);

    std::vector<long long> prime_full;
    prime_full.reserve(n-1);

    std::vector<bool> marks(limit + 1, false);

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
            prime.push_back(p);
        }
    }


    long long low = limit;
    long long high = 2*limit;


    while (low < n)
    {
        if (high >= n)
           high = n;


        bool mark[limit+1];
        memset(mark, false, sizeof(mark));


        for (long long i = 0; i < prime.size(); i++)
        {

            long long loLim = floor(low/prime[i]) * prime[i];
            if (loLim < low)
                loLim += prime[i];
 

            for (long long j=loLim; j<high; j+=prime[i])
                mark[j-low] = true;
        }

        for (long long i = low; i<high; i++)
            if (mark[i - low] == false)
                prime_full.push_back(i);


        low = low + limit;
        high = high + limit;
    }


    auto t2 = std::chrono::high_resolution_clock::now();

    // Output Execution time
    outputFile << "Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0 << "s" << std::endl
               << std::endl;
    std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0 << "s" << std::endl
              << std::endl;

    for (auto v:prime){
        outputFile << v<< std::endl;
    }
    for (auto v:prime_full){
        outputFile << v<< std::endl;
    }


}

int main(int argc, char *argv[])
{
    long long n;

    std::cout << "N: ";
    std::cin >> n;  
    std::cout << std::endl;

    // std::ofstream outputFile1;
    // outputFile1.open("sieveoferatosthenes.txt");
    // SieveOfEratosthenes(n, outputFile1);
    // outputFile1.close();

    std::ofstream outputFile2;
    outputFile2.open("sieveoferatosthenesfastmarking.txt");
    SieveOfEratosthenesFastMarking(n, outputFile2);
    outputFile2.close();

    std::ofstream outputFile3;
    outputFile3.open("sieveoferatosthenesfastmarkingreorganized.txt");
    SieveOfEratosthenesFastMarkingReorganized(n, outputFile3);
    outputFile3.close();
}