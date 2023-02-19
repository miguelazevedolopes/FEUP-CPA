#include <iostream>
#include <fstream>
#include <eml.h>
#include <omp.h>

using namespace std;

#define SYSTEMTIME clock_t

void OnMult(int m_ar, int m_br, ofstream &out)
{

    SYSTEMTIME Time1, Time2;


    char st[200];
    double temp;
    int i, j, k;

    out << m_ar<<"x"<<m_br<<"Simple:"<< endl;

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

    size_t count;
    emlDeviceGetCount(&count);
    emlData_t* data[count];

    Time1 = clock();
    emlStart();

    #pragma omp parallel for
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
    emlStop(data);
    Time2 = clock();

    double consumed, elapsed;
    emlDataGetConsumed(data[0], &consumed);
    emlDataGetElapsed(data[0], &elapsed);
    emlDataFree(data[0]);

    sprintf(st,"This device consumed %g J in %g s\n", consumed, elapsed);
    out << st;

    sprintf(st, "Time: %3.3f seconds\n", (double)(Time2 - Time1) / CLOCKS_PER_SEC);
    out << st;

    // display 10 elements of the result matrix tto verify correctness
    out << "Result matrix: " << endl;
    for (i = 0; i < 1; i++)
    {
        for (j = 0; j < min(10, m_br); j++)
            out << phc[j] << " ";
    }
    out << endl<<endl;
    cout<<"End mult matrix"<<endl;

    free(pha);
    free(phb);
    free(phc);
}

// add code here for line x line matriz multiplication
void OnMultLine(int m_ar, int m_br, ofstream &out)
{
    SYSTEMTIME Time1, Time2;

    char st[200];
    double temp;
    int i, j, k;

    double *pha, *phb, *phc;

    pha = (double *)malloc((m_ar * m_ar) * sizeof(double));
    phb = (double *)malloc((m_ar * m_ar) * sizeof(double));
    phc = (double *)malloc((m_ar * m_ar) * sizeof(double));

    out << m_ar<<"x"<<m_br<<"MultLine:"<< endl;


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

    size_t count;
    emlDeviceGetCount(&count);
    emlData_t* data[count];

    Time1 = clock();
    emlStart();
    #pragma omp parallel for
    for (i = 0; i < m_ar; i++) {
        for (k = 0; k < m_br; k++) {
            for (j = 0; j < m_ar; j++) {
                phc[i*m_ar + j] += pha[i*m_ar + k] * phb[k*m_br + j];
            }
        }
    }

    emlStop(data);
    Time2 = clock();

    double consumed, elapsed;
    emlDataGetConsumed(data[0], &consumed);
    emlDataGetElapsed(data[0], &elapsed);
    emlDataFree(data[0]);

    sprintf(st,"This device consumed %g J in %g s\n", consumed, elapsed);
    out << st;

    sprintf(st, "Time: %3.3f seconds\n", (double)(Time2 - Time1) / CLOCKS_PER_SEC);
    out << st;

    // display 10 elements of the result matrix tto verify correctness
    out << "Result matrix: " << endl;
    for (i = 0; i < 1; i++)
    {
        for (j = 0; j < min(10, m_br); j++)
            out << phc[j] << " ";
    }
    out << endl<<endl;
    cout<<"End multline matrix"<<endl;

    free(pha);
    free(phb);
    free(phc);

}

int main(int argc, char *argv[])
{
    emlInit();
    ofstream outputFile;
    outputFile.open("output2.txt");


    OnMult(1000, 1000, outputFile);

    OnMult(1400, 1400, outputFile);

    OnMult(2000, 2000, outputFile);

    OnMultLine(1000, 1000, outputFile);

    OnMultLine(1400, 1400, outputFile);

    OnMultLine(2000, 2000, outputFile);

    outputFile.close();
    emlShutdown();
}