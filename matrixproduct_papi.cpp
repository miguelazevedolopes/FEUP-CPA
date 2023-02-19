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

    
    printf("1000x1000 Simple:\n");
    ret = PAPI_start(EventSet);
    if (ret != PAPI_OK)
        printf("ERROR: Start PAPI\n");
    OnMult(1000, 1000,outputFile);
    ret = PAPI_stop(EventSet, values);
    if (ret != PAPI_OK)
        printf("ERROR: Stop PAPI\n");

    printf("PAPI_L1_DCM: %llu\n",values[0]);
    printf("PAPI_LD_INS: %llu\n",values[0]);
    printf("PAPI_SR_INS: %llu\n\n",values[0]);


    printf("1400x1400 Simple:\n");
    ret = PAPI_start(EventSet);
    if (ret != PAPI_OK)
        printf("ERROR: Start PAPI\n");
    OnMult(1400, 1400,outputFile);
    ret = PAPI_stop(EventSet, values);
    if (ret != PAPI_OK)
        printf("ERROR: Stop PAPI\n");

    printf("PAPI_L1_DCM: %llu\n",values[0]);
    printf("PAPI_LD_INS: %llu\n",values[0]);
    printf("PAPI_SR_INS: %llu\n\n",values[0]);


    printf("2000x2000 Simple:\n");
    ret = PAPI_start(EventSet);
    if (ret != PAPI_OK)
        printf("ERROR: Start PAPI\n");
    OnMult(2000, 2000,outputFile);
    ret = PAPI_stop(EventSet, values);
    if (ret != PAPI_OK)
        printf("ERROR: Stop PAPI\n");

    printf("PAPI_L1_DCM: %llu\n",values[0]);
    printf("PAPI_LD_INS: %llu\n",values[0]);
    printf("PAPI_SR_INS: %llu\n\n",values[0]);


    printf("1000x1000 MultLine:\n");
    ret = PAPI_start(EventSet);
    if (ret != PAPI_OK)
        printf("ERROR: Start PAPI\n");
    OnMultLine(1000, 1000,outputFile);
    ret = PAPI_stop(EventSet, values);
    if (ret != PAPI_OK)
        printf("ERROR: Stop PAPI\n");

    printf("PAPI_L1_DCM: %llu\n",values[0]);
    printf("PAPI_LD_INS: %llu\n",values[0]);
    printf("PAPI_SR_INS: %llu\n\n",values[0]);


    printf("1400x1400 MultLine:\n");
    ret = PAPI_start(EventSet);
    if (ret != PAPI_OK)
        printf("ERROR: Start PAPI\n");
    OnMultLine(1400, 1400,outputFile);
    ret = PAPI_stop(EventSet, values);
    if (ret != PAPI_OK)
        printf("ERROR: Stop PAPI\n");

    printf("PAPI_L1_DCM: %llu\n",values[0]);
    printf("PAPI_LD_INS: %llu\n",values[0]);
    printf("PAPI_SR_INS: %llu\n\n",values[0]);


    printf("2000x2000 MultLine:\n");
    ret = PAPI_start(EventSet);
    if (ret != PAPI_OK)
        printf("ERROR: Start PAPI\n");
    OnMultLine(2000, 2000,outputFile);
    ret = PAPI_stop(EventSet, values);
    if (ret != PAPI_OK)
        printf("ERROR: Stop PAPI\n");

    printf("PAPI_L1_DCM: %llu\n",values[0]);
    printf("PAPI_LD_INS: %llu\n",values[0]);
    printf("PAPI_SR_INS: %llu\n\n",values[0]);


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

    outputFile.close();
}