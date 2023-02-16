#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <cstdlib>
#include <papi.h>
#include <fstream>

using namespace std;

#define SYSTEMTIME clock_t

void OnMult(int m_ar, int m_br, ofstream &out)
{

    SYSTEMTIME Time1, Time2;

    char st[100];
    double temp;
    int i, j, k;

    double *pha, *phb, *phc;

    pha = (double *)malloc((m_ar * m_ar) * sizeof(double));
    phb = (double *)malloc((m_ar * m_ar) * sizeof(double));
    phc = (double *)malloc((m_ar * m_ar) * sizeof(double));

    for (i = 0; i < m_ar; i++)
        for (j = 0; j < m_ar; j++)
            pha[i * m_ar + j] = (double)1.0;

    for (i = 0; i < m_br; i++)
        for (j = 0; j < m_br; j++)
            phb[i * m_br + j] = (double)(i + 1);

    Time1 = clock();

    for (i = 0; i < m_ar; i++)
    {
        for (j = 0; j < m_br; j++)
        {
            temp = 0;
            for (k = 0; k < m_ar; k++)
            {
                temp += pha[i * m_ar + k] * phb[k * m_br + j];
            }
            phc[i * m_ar + j] = temp;
        }
    }

    Time2 = clock();
    sprintf(st, "Time: %3.3f seconds\n", (double)(Time2 - Time1) / CLOCKS_PER_SEC);
    out << st;

    // display 10 elements of the result matrix tto verify correctness
    out << "Result matrix: " << endl;
    for (i = 0; i < 1; i++)
    {
        for (j = 0; j < min(10, m_br); j++)
            out << phc[j] << " ";
    }
    out << endl;
    cout<<"End mult matrix"<<endl;

    free(pha);
    free(phb);
    free(phc);
}

// add code here for line x line matriz multiplication
void OnMultLine(int m_ar, int m_br, ofstream &out)
{
    SYSTEMTIME Time1, Time2;

    char st[100];
    double temp;
    int i, j, k;

    double *pha, *phb, *phc;

    pha = (double *)malloc((m_ar * m_ar) * sizeof(double));
    phb = (double *)malloc((m_ar * m_ar) * sizeof(double));
    phc = (double *)malloc((m_ar * m_ar) * sizeof(double));

    for (i = 0; i < m_ar; i++){
        for (j = 0; j < m_ar; j++){
            pha[i * m_ar + j] = (double)1.0;
            phc[i * m_ar + j] = (double)0.0;
        }
    }

    for (i = 0; i < m_br; i++){
        for (j = 0; j < m_br; j++){
            phb[i * m_br + j] = (double)(i + 1);
        }
    }

    Time1 = clock();

    for (i = 0; i < m_ar; i++) {
        for (k = 0; k < m_br; k++) {
            for (j = 0; j < m_ar; j++) {
                phc[i*m_ar + j] += pha[i*m_ar + k] * phb[k*m_br + j];
            }
        }
    }

    Time2 = clock();
    sprintf(st, "Time: %3.3f seconds\n", (double)(Time2 - Time1) / CLOCKS_PER_SEC);
    out << st;

    // display 10 elements of the result matrix tto verify correctness
    out << "Result matrix: " << endl;
    for (i = 0; i < 1; i++)
    {
        for (j = 0; j < min(10, m_br); j++)
            out << phc[j] << " ";
    }
    out << endl;
    cout<<"End multline matrix"<<endl;

    free(pha);
    free(phb);
    free(phc);

}

// add code here for block x block matriz multiplication
void OnMultBlock(int m_ar, int m_br, int bkSize, ofstream &out)
{
    SYSTEMTIME Time1, Time2;

    char st[100];
    double temp;
    int i, j, k, jj, kk;

    double *pha, *phb, *phc;

    pha = (double *)malloc((m_ar * m_ar) * sizeof(double));
    phb = (double *)malloc((m_ar * m_ar) * sizeof(double));
    phc = (double *)malloc((m_ar * m_ar) * sizeof(double));

    for (i = 0; i < m_ar; i++){
        for (j = 0; j < m_ar; j++){
            pha[i * m_ar + j] = (double)1.0;
            phc[i * m_ar + j] = (double)0.0;
        }
    }

    for (i = 0; i < m_br; i++){
        for (j = 0; j < m_br; j++){
            phb[i * m_br + j] = (double)(i + 1);
        }
    }

    Time1 = clock();

    for (kk = 0; kk < m_ar; kk += bkSize){
        for (jj = 0; jj < m_ar; jj += bkSize){
            for (i = 0 ; i < m_ar; i++){
                for (k = kk; k < kk + bkSize; k++){
                    for (j = jj ; j < jj + bkSize; j++){
                       phc[i*m_ar + j] += pha[i*m_ar + k] * phb[k*m_br + j] ;
                    }
                }
            }
        }
    }

    Time2 = clock();
    sprintf(st, "Time: %3.3f seconds\n", (double)(Time2 - Time1) / CLOCKS_PER_SEC);
    out << st;

    // display 10 elements of the result matrix tto verify correctness
    out << "Result matrix: " << endl;
    for (i = 0; i < 1; i++)
    {
        for (j = 0; j < min(10, m_br); j++)
            out << phc[j] << " ";
    }
    out << endl;

    free(pha);
    free(phb);
    free(phc);

    cout<<"End block matrix"<<endl;
}

void handle_error(int retval)
{
    printf("PAPI error %d: %s\n", retval, PAPI_strerror(retval));
    exit(1);
}

void init_papi()
{
    int retval = PAPI_library_init(PAPI_VER_CURRENT);
    if (retval != PAPI_VER_CURRENT && retval < 0)
    {
        printf("PAPI library version mismatch!\n");
        exit(1);
    }
    if (retval < 0)
        handle_error(retval);

    std::cout << "PAPI Version Number: MAJOR: " << PAPI_VERSION_MAJOR(retval)
              << " MINOR: " << PAPI_VERSION_MINOR(retval)
              << " REVISION: " << PAPI_VERSION_REVISION(retval) << "\n";
}

int main(int argc, char *argv[])
{
    char c;
    int op;
    ofstream outputFile;
    outputFile.open("output1.txt");
    int EventSet = PAPI_NULL;
    long long values[3];
    int ret;

    ret = PAPI_library_init(PAPI_VER_CURRENT);
    if (ret != PAPI_VER_CURRENT)
        outputFile << "FAIL" << endl;

    ret = PAPI_create_eventset(&EventSet);
    if (ret != PAPI_OK)
        outputFile << "ERROR: create eventset" << endl;

    ret = PAPI_add_event(EventSet, PAPI_L1_DCM);
    if (ret != PAPI_OK)
        outputFile << "ERROR: PAPI_L1_DCM" << endl;

    ret = PAPI_add_event(EventSet, PAPI_LD_INS);
    if (ret != PAPI_OK)
        outputFile << "ERROR: PAPI_LD_INS" << endl;

    ret = PAPI_add_event(EventSet, PAPI_SR_INS);
    if (ret != PAPI_OK){
        outputFile << "ERROR: PAPI_SR_INS" << endl;
    }


    for(int i=600;i<=3000;i+=400){
        ret = PAPI_start(EventSet);
        if (ret != PAPI_OK)
            outputFile << "ERROR: Start PAPI" << endl;
        outputFile << "Matrix size: "<<i<<endl;
        OnMult(i, i,outputFile);
        ret = PAPI_stop(EventSet, values);
        if (ret != PAPI_OK)
            outputFile << "ERROR: Stop PAPI" << endl;

        outputFile << values[0]<<endl;
        outputFile <<values[1]<<endl;
        outputFile <<values[2]<<endl;
        outputFile <<"--------"<<endl;
    }

    

    for(int i=600;i<=3000;i+=400){
        ret = PAPI_start(EventSet);
        if (ret != PAPI_OK)
            outputFile << "ERROR: Start PAPI" << endl;
        outputFile << "Matrix size: "<<i<<endl;
        OnMultLine(i, i,outputFile);
        ret = PAPI_stop(EventSet, values);
        if (ret != PAPI_OK)
            outputFile << "ERROR: Stop PAPI" << endl;

        outputFile <<values[0]<<endl;
        outputFile << values[1]<<endl;
        outputFile <<values[2]<<endl;
        outputFile <<"--------"<<endl;
    }

    

    
    for(int i=4096;i<=10240;i+=2048){
        ret = PAPI_start(EventSet);
        if (ret != PAPI_OK)
            outputFile << "ERROR: Start PAPI" << endl;

        outputFile << "Matrix size: "<<i<<endl;
        OnMultLine(i, i,outputFile);
        ret = PAPI_stop(EventSet, values);
        if (ret != PAPI_OK)
            outputFile << "ERROR: Stop PAPI" << endl;

        outputFile << values[0]<<endl;
        outputFile << values[1]<<endl;
        outputFile << values[2]<<endl;
        outputFile <<"--------"<<endl;

    }

    
   
    for(int i=4096;i<=10240;i+=2048){
        ret = PAPI_start(EventSet);
        if (ret != PAPI_OK)
            outputFile << "ERROR: Start PAPI" << endl;
        outputFile << "Matrix size: "<<i<<endl;
        outputFile<< "Block size:"<< 128<<endl;
        OnMultBlock(i, i, 128,outputFile);
        ret = PAPI_stop(EventSet, values);
        if (ret != PAPI_OK)
            outputFile << "ERROR: Stop PAPI" << endl;

        outputFile << values[0]<<endl;
        outputFile << values[1]<<endl;
        outputFile << values[2]<<endl;
        outputFile <<"--------"<<endl;
    }

    
    
    

    for(int i=4096;i<=10240;i+=2048){
        ret = PAPI_start(EventSet);
        if (ret != PAPI_OK)
            outputFile << "ERROR: Start PAPI" << endl;
        outputFile << "Matrix size: "<<i<<endl;
        outputFile<< "Block size:"<< 256<<endl;
        OnMultBlock(i, i, 256,outputFile);
        ret = PAPI_stop(EventSet, values);
        if (ret != PAPI_OK)
            outputFile << "ERROR: Stop PAPI" << endl;

        outputFile << values[0]<<endl;
        outputFile << values[1]<<endl;
        outputFile << values[2]<<endl;
        outputFile <<"--------"<<endl;
    }

    


    for(int i=4096;i<=10240;i+=2048){
        ret = PAPI_start(EventSet);
        if (ret != PAPI_OK)
            outputFile << "ERROR: Start PAPI" << endl;

        outputFile << "Matrix size: "<<i<<endl;
        outputFile<< "Block size:"<< 512<<endl;
        OnMultBlock(i, i, 512,outputFile);
        ret = PAPI_stop(EventSet, values);
        if (ret != PAPI_OK)
            outputFile << "ERROR: Stop PAPI" << endl;

        outputFile <<  values[0]<<endl;
        outputFile << values[1]<<endl;
        outputFile << values[2]<<endl;
        outputFile <<"--------"<<endl;

    }

    
    
    for(int i=4096;i<=10240;i+=2048){
        ret = PAPI_start(EventSet);
        if (ret != PAPI_OK)
            outputFile << "ERROR: Start PAPI" << endl;
        outputFile << "Matrix size: "<<i<<endl;
        outputFile<< "Block size:"<< 1024<<endl;
        OnMultBlock(i, i, 1024,outputFile);
        ret = PAPI_stop(EventSet, values);
        if (ret != PAPI_OK)
            outputFile << "ERROR: Stop PAPI" << endl;

        outputFile << values[0]<<endl;
        outputFile << values[1]<<endl;
        outputFile << values[2]<<endl;
        outputFile <<"--------"<<endl;
    
    }

    
    ret = PAPI_reset(EventSet);
    if (ret != PAPI_OK)
        outputFile << "FAIL reset" << endl;

    ret = PAPI_remove_event(EventSet, PAPI_L1_DCM);
    if (ret != PAPI_OK)
        outputFile << "FAIL remove event" << endl;

    ret = PAPI_remove_event(EventSet, PAPI_L2_DCM);
    if (ret != PAPI_OK)
        outputFile << "FAIL remove event" << endl;

    ret = PAPI_remove_event(EventSet, PAPI_LST_INS);
    if (ret != PAPI_OK)
        outputFile << "FAIL remove event" << endl;
    

    ret = PAPI_destroy_eventset(&EventSet);
    if (ret != PAPI_OK)
        outputFile << "FAIL destroy" << endl;
}