#include <iostream>
#include <ctime>
#include <algorithm>
#include <cmath>
#include <fstream>


using namespace std;
const double pi = 3.141592653589793;

int part(float p[], int in[], int start, int _end)
{
    float pivot = p[_end];
    int pIndex = start;
    for (int i = start; i < _end; i++)
    {
        if (p[i] <= pivot)
        {
            swap(p[i], p[pIndex]);
            swap(in[i], in[pIndex]);
            pIndex++;
        }
    }
    swap (p[pIndex], p[_end]);
    swap (in[pIndex], in[_end]);
    return pIndex;
}
void quicksort(float p[], int in[], int start, int _end)
{
    if (start >= _end) {
        return;
    }
    int pivot = part(p, in, start, _end);
    quicksort(p, in, start, pivot - 1);
    quicksort(p, in, pivot + 1, _end);
}
int part(float p[], float in[], int start, int _end)
{
    float pivot = p[_end];
    int pIndex = start;
    for (int i = start; i < _end; i++)
    {
        if (p[i] <= pivot)
        {
            swap(p[i], p[pIndex]);
            swap(in[i], in[pIndex]);
            pIndex++;
        }
    }
    swap (p[pIndex], p[_end]);
    swap (in[pIndex], in[_end]);
    return pIndex;
}
void quicksort(float p[], float in[], int start, int _end)
{
    if (start >= _end) {
        return;
    }
    int pivot = part(p, in, start, _end);
    quicksort(p, in, start, pivot - 1);
    quicksort(p, in, pivot + 1, _end);
}
int part(int p[], float in[], int start, int _end)
{
    int pivot = p[_end];
    int pIndex = start;
    for (int i = start; i < _end; i++)
    {
        if (p[i] <= pivot)
        {
            swap(p[i], p[pIndex]);
            swap(in[i], in[pIndex]);
            pIndex++;
        }
    }
    swap (p[pIndex], p[_end]);
    swap (in[pIndex], in[_end]);
    return pIndex;
}
void quicksort(int p[], float in[], int start, int _end)
{
    if (start >= _end) {
        return;
    }
    int pivot = part(p, in, start, _end);
    quicksort(p, in, start, pivot - 1);
    quicksort(p, in, pivot + 1, _end);
}

void cross(int fights, int m, int cross_switch, int **children, int **parents)
{
    int i, j, k, p, splitter, splitter2;
    if (cross_switch == 1)
        {
            for (i = 0, k = 0; k < fights; i++, k += 2)
            {
                splitter = 1 + (rand() % (m-1));
                for (j = 0; j < splitter; j++)
                {
                    children[i][j] = parents[k][j];
                }
                for (j = splitter; j < m; j++)
                {
                    children[i][j] = parents[k + 1][j];
                }
            }
        }
        if (cross_switch == 2)
        {
            for (i = 0, k = 0; k < fights; i++, k+=2)
            {
                splitter = 1 + (rand() % (m-1));
                splitter2 = 1 + (rand() % (m-1));
                if (splitter < splitter2)
                {
                     for (j = 0; j < splitter; j++)
                    {
                        children[i][j] = parents[k][j];
                    }
                    for (j = splitter; j < splitter2; j++)
                    {
                        children[i][j] = parents[k + 1][j];
                    }
                    for (j = splitter2; j < m; j++)
                    {
                        children[i][j] = parents[k][j];
                    }
                }
                if (splitter2 < splitter)
                {
                     for (j = 0; j < splitter2; j++)
                    {
                        children[i][j] = parents[k][j];
                    }
                    for (j = splitter2; j < splitter; j++)
                    {
                        children[i][j] = parents[k + 1][j];
                    }
                    for (j = splitter; j < m; j++)
                    {
                        children[i][j] = parents[k][j];
                    }
                }
               if (splitter == splitter2)
               {
                   k-=2;
                   i--;
                   continue;
               }
            }
        }
        if (cross_switch == 3)
        {
            for (i = 0, k = 0; k < fights; i++, k+=2)
            {
                for (j = 0; j < m; j++)
                {
                    p = rand()%2;
                    children[i][j] = parents[k+p][j];
                }
            }
        }
}
void mutation(string mut_switch, int n, int m, int **children)
{
    float mut, p_weak = 1./(3.*m), p_average = 1./m, p_strong = 3./m;
    int i, j;
    if(mut_switch == "weak"){

            for (i = 0;i< n;i++)
            {
                for (j = 0; j<m;j++)
                {
                     mut = ((double) rand() / (double)(RAND_MAX));
                    if (mut <= p_weak)
                    {
                        if (children[i][j] == 0)
                            children[i][j] = 1;
                        if (children[i][j] == 1)
                            children[i][j] = 0;
                    }
                }
            }
        }
        if(mut_switch == "average"){

            for (i = 0;i< n;i++)
            {
                for (j = 0; j<m;j++)
                {
                     mut = ((double) rand() / (double)(RAND_MAX));
                    if (mut <= p_average)
                    {
                        if (children[i][j] == 0)
                            children[i][j] = 1;
                        if (children[i][j] == 1)
                            children[i][j] = 0;
                    }
                }
            }
        }
        if(mut_switch == "strong")
        {
            for (i = 0;i< n;i++)
            {
                for (j = 0; j<m;j++)
                {
                     mut = ((double) rand() / double(RAND_MAX));
                    if (mut <= p_strong)
                    {
                        if (children[i][j] == 0)
                            children[i][j] = 1;
                        if (children[i][j] == 1)
                            children[i][j] = 0;
                    }
                }
            }
        }
}

void selection (string sel_switch, float *power, int fights, int n, int **parents, int **A, int m)
{
    float *temp_power = new float[n];
    float sumpower = 0;
    int *marks = new int[n];
    int *index = new int[n];
    int fight, i, j, candidate1, candidate2, winner, cont = 0, k, sumrang = 0;
    float *rang = new float[n];
    float *p = new float[n];
    float temp = 0, randp;
    if (sel_switch == "prop")
    {
        for (i = 0; i < n; i++)
        {
            sumpower += power[i];
        }
        for (i = 0; i < n; i++)
        {
            p[i] = temp + (float)power[i]/((float)sumpower);
            temp += (float)power[i]/((float)sumpower);
        }
        for (fight = 0; fight < fights; fight++)
        {
            randp = ((double) rand() / (RAND_MAX));
            if (randp >= 0&&randp < p[0])
            {
                for (i = 0; i < m; i++)
                    parents[fight][i] = A[0][i];
                continue;
            }
            for (j = 0; j < n-1; j++)
            {
                if (randp < p[j+1]&&randp >=p[j])
                {
                    for (i = 0; i < m; i++)
                        parents[fight][i] = A[j+1][i];
                    break;
                }
            }
        }

    }
    if (sel_switch == "rang")
    {
        for (i = 0; i < n; i++)
        {
            temp_power[i] = power[i];
            index[i] = i;
            rang[i] = i + 1;
            marks[i] = 0;
        }
        quicksort(temp_power, index, 0, n - 1);
        //если пригодности одинаковые
        for (i = 0; i < n-1; i++)
        {
            if (marks[i] == 2)
                continue;
            sumrang = i+1;
            for (j = i + 1; j < n; j++)
            {
                if (temp_power[i] == temp_power[j])
                {
                    cont +=1;
                    marks[i] = 1;
                    marks[j] = 1;
                    sumrang += j+1;
                }
            }
            if (cont > 0)
            {
                cont += 1;
                for (k = 0; k < n; k++)
                {
                    if (marks[k] == 1)
                    {
                        rang[k] = (float)sumrang/((float) cont);
                        marks[k] = 2;
                    }
                }
            }
            cont = 0;
            sumrang = 0;
        }

        quicksort(index, rang, 0, n-1);

        for (i = 0; i < n; i++)
        {
            sumpower += rang[i];
        }
        for (i = 0; i < n; i++)
        {
            p[i] = temp + (float)rang[i]/((float)sumpower);
            //cout << p[i] << endl;
            temp += (float)rang[i]/((float)sumpower);
        }

        for (fight = 0; fight < fights; fight++)
        {
            randp = ((double) rand() / (RAND_MAX));
            if (0 <= randp < p[0])
            {
                for (i = 0; i < m; i++)
                    parents[fight][i] = A[0][i];
                continue;
            }
            for (j = 0; j < n-1; j++)
            {
                if (p[j] <= randp < p[j+1])
                {
                    for (i = 0; i < m; i++)
                        parents[fight][i] = A[j+1][i];
                    break;
                }
            }
        }

    }
    if (sel_switch == "tour")
    {
         for (int fight = 0; fight < fights; fight++)
        {
            candidate1 = rand() % n;
            candidate2 = rand() % n;
            if (power[candidate1] > power[candidate2])
            {
                for (i = 0; i < m; i++)
                    parents[fight][i] = A[candidate1][i];
            }
            if (power[candidate1] < power[candidate2])
            {
                for (i = 0; i < m; i++)
                    parents[fight][i] = A[candidate2][i];
            }
            if (power[candidate1] == power[candidate2])
            {
                winner = rand() % 2;
                if (winner == 0)
                    for (i = 0; i < m; i++)
                        parents[fight][i] = A[candidate1][i];
                if (winner == 1)
                    for (i = 0; i < m; i++)
                        parents[fight][i] = A[candidate2][i];

            }
        }
    }
    delete[] marks;
    delete[] temp_power;
    delete[] index;
    delete[] rang;
    delete[] p;
}



void count_x(int **A, int n, float *h, float *left, float **x, int *number_m, int number_x, bool parents = false)
{
    int temp = 0, i, j, k, t, cum_m=0;
    for (i = 0; i < n; i++)
    {
        cum_m = 0;
        for (t = 0; t < number_x; t++)
        {
            cum_m += number_m[t];
            if (t == 0)
            {
                for (j = 0, k = number_m[t]-1; j < number_m[t]; j++, k--)
                {
                    if (A[i][j] == 1)
                    {
                        temp += pow(2,k);
                    }
                }
            }
            else
            {
                for (j = cum_m-number_m[t], k = number_m[t]-1; j < cum_m-1; j++, k--)
                {
                    if (A[i][j] == 1)
                    {
                        temp += pow(2,k);
                    }
                }
            }
            if (parents)
                x[i][t] = h[t]*temp + left[t];
            else
                x[i+n][t] = h[t]*temp + left[t];
            temp = 0;
        }

    }
}

void func(float **x, float *f, int n, int kind_of_func, int number_x, bool parents = false)
{
    int i, j = 0, k;
    float sum1, sum2;
    if (parents == false)
    {
        j = n;
        n = n*2;
    }
    if (kind_of_func == 0)
        for (i = j; i < n; i++)
            f[i] = 10*x[i][0]-5*x[i][1];
    if (kind_of_func == 1)
        for (i = j; i < n; i++)
            f[i] = x[i][0]*x[i][0]+x[i][1]*x[i][1];
    if (kind_of_func == 2)
        for (i = j; i < n; i++)
            f[i] = 5*x[i][0]+0.5*x[i][1];
    if (kind_of_func == 3)
        for (i = j; i < n; i++)
            f[i] = x[i][0]*x[i][0]*x[i][0]+x[i][1}*x[i][1]+x[i][2];
    if (kind_of_func == 4)
        for (i = j; i < n; i++)
            f[i] = 2000*x[i][0]+2400*x[i][1];
    if (kind_of_func == 5)
    {
        for (i = j; i < n; i++)
        {
            sum1 = 0;
            sum2 = 0;
            for (k = 0; k < number_x; k++)
            {
                sum1+=x[i][k]*x[i][k]/number_x;
                sum2+=cos(2*pi*x[i][k]/number_x);
            }
            // e это эпсилант??
            f[i] = 20-20*exp(-0.2*sum1)-exp(sum2);
        }
    }
    if (kind_of_func == 6)
    {
        for (i = j; i < n; i++)
        {
            sum1 = 0;
            for (k = 0; k < number_x; k++)
                sum1+=0.1*x[i][k]*x[i][k]-4*cos(0.8*x[i][k])+4;
            // e это эпсилант??
            f[i] = sum1;
        }
    }
}

void fine(float **x, int n, float *fines, float b, int kind_of_func, bool parents = false)
{
    int i, j, k = 0;
    float mult_fines = 1, temp_fines = 0;
    if (parents == false)
    {
        k = n;
        n = n*2;
    }
    if(kind_of_func == 0)
    {
        for(i = k; i < n; i++)
        {
            temp_fines = 0;
            for (j = 0; j < b; j++)
                mult_fines*=max(float(0), x[i][1]-15);
            temp_fines+=mult_fines;
            mult_fines = 1;
            for (j = 0; j < b; j++)
                mult_fines*=max(float(0), x[i][1]+2*x[i][0]*x[i][0]-20);
            temp_fines+=mult_fines;
            mult_fines = 1;
            for (j = 0; j < b; j++)
                mult_fines*=max(0., -((x[i][0]*x[i][0])/2.)-x[i][1]);
            temp_fines+=mult_fines;
            fines[i] = temp_fines;
        }

    }
    if(kind_of_func == 1)
    {
        for(i = k; i < n; i++)
        {
            temp_fines = 0;
            for (j = 0; j < b; j++)
                mult_fines*=max(float(0), x[i][1]-7-sin(2*x[i][0]));
            temp_fines+=mult_fines;
            mult_fines = 1;
            for (j = 0; j < b; j++)
                mult_fines*=max(float(0), 1-sin(2*x[i][0])-x[i][1]);
            temp_fines+=mult_fines;
            mult_fines = 1;
            if (x[i][0]>=0 && x[i][0]<=4)
                mult_fines = 0;
            if (x[i][0]<0)
                for (j = 0; j < b; j++)
                    mult_fines*=-x[i][0];
            if (x[i][0]>4)
                for (j = 0; j < b; j++)
                    mult_fines*=x[i][0]-4;
            temp_fines+=mult_fines;
            fines[i] = temp_fines;
        }

    }
    if(kind_of_func == 2)
    {
        for(i = k; i < n; i++)
        {
            temp_fines = 0;
            for (j = 0; j < b; j++)
                mult_fines*=max(float(0), x[i][1]+2*x[i][0]-5);
            temp_fines+=mult_fines;
            mult_fines = 1;
            for (j = 0; j < b; j++)
                mult_fines*=max(0., x[i][0]-1.5-x[i][1]);
            temp_fines+=mult_fines;
            mult_fines = 1;
            for (j = 0; j < b; j++)
                mult_fines*=max(float(0), x[i][1]-2*x[i][0]-1);
            temp_fines+=mult_fines;
            fines[i] = temp_fines;
        }

    }
    if(kind_of_func == 3)
    {
        for(i = k; i < n; i++)
        {
            temp_fines = 0;
            for (j = 0; j < b; j++)
                mult_fines*=max(float(0), -x[i][0]);
            temp_fines+=mult_fines;
            mult_fines = 1;
            for (j = 0; j < b; j++)
                mult_fines*=max(0., -x[i][1]);
            temp_fines+=mult_fines;
            mult_fines = 1;
            for (j = 0; j < b; j++)
                mult_fines*=max(float(0), -x[i][2]);
            temp_fines+=mult_fines;
            mult_fines = 1;
            for (j = 0; j < b; j++)
                mult_fines*=max(float(0), x[i][0]*x[i][0]+x[i][1]*x[i][1]+x[i][2]*x[i][2]-25);
            temp_fines+=mult_fines;
            fines[i] = temp_fines;
        }

    }
    if(kind_of_func == 4)
    {
        for(i = k; i < n; i++)
        {
            temp_fines = 0;
            for (j = 0; j < b; j++)
                mult_fines*=max(float(0), -x[i][0]);
            temp_fines+=mult_fines;
            mult_fines = 1;
            for (j = 0; j < b; j++)
                mult_fines*=max(0., -x[i][1]);
            temp_fines+=mult_fines;
            mult_fines = 1;
            for (j = 0; j < b; j++)
                mult_fines*=max(float(0), x[i][0]/120+x[i][1]/110-1);
            temp_fines+=mult_fines;
            mult_fines = 1;
            for (j = 0; j < b; j++)
                mult_fines*=max(float(0), 4*x[i][0]+x[i][1]-320);
            temp_fines+=mult_fines;
            mult_fines = 1;
            for (j = 0; j < b; j++)
                mult_fines*=max(float(0), x[i][0]+x[i][1]-110);
            temp_fines+=mult_fines;
            mult_fines = 1;
            for (j = 0; j < b; j++)
                mult_fines*=max(float(0), x[i][0]/340+x[i][1]/120-1);
            temp_fines+=mult_fines;
            mult_fines = 1;
            for (j = 0; j < b; j++)
                mult_fines*=max(float(0), x[i][0]+2*x[i][1]-160);
            temp_fines+=mult_fines;
            mult_fines = 1;
            for (j = 0; j < b; j++)
                mult_fines*=max(float(0), x[i][0]+4*x[i][1]-280);
            temp_fines+=mult_fines;
            fines[i] = temp_fines;
        }

    }
    if(kind_of_func == 5)
    {
        for(i = k; i < n; i++)
        {
            temp_fines = 0;
            for (j = 0; j < b; j++)
                mult_fines*=max(float(0), 2*x[i][0]-3*x[i][1]+4*x[i][2]-10);
            temp_fines+=mult_fines;
            mult_fines = 1;
            for (j = 0; j < b; j++)
                mult_fines*=max(0., 4*x[i][1]-5*x[i][2]+x[i][3]-1);
            temp_fines+=mult_fines;
            mult_fines = 1;
            for (j = 0; j < b; j++)
                mult_fines*=max(float(0), 10*x[i][0]+7.5*x[i][2]-8.4*x[i][3]-3.5);
            temp_fines+=mult_fines;
            mult_fines = 1;
            for (j = 0; j < b; j++)
                mult_fines*=max(float(0), -3.1*x[i][0]+21.7*x[i][1]-36.4*x[i][3]-16.2);
            temp_fines+=mult_fines;
            fines[i] = temp_fines;
        }

    }
    if(kind_of_func == 6)
    {
        for(i = k; i < n; i++)
        {
            temp_fines = 0;
            for (j = 0; j < b; j++)
                mult_fines*=max(float(0), x[i][0]*x[i][0]+9*x[i][1]*x[i][1]-36);
            temp_fines+=mult_fines;
            mult_fines = 1;
            for (j = 0; j < b; j++)
                mult_fines*=max(float(0), 9*x[i][0]*x[i][0]+x[i][1]*x[i][1]-36);
            temp_fines+=mult_fines;
            fines[i] = temp_fines;
        }

    }
}
//нужно сохранять значение лямбды предыдущее
void fitness(float *power, float* f, int globali, float *fines, float C, float a, float b, int n, int kind_of_func,
             int kind_of_fine, float b1, float b2, float k, int without_fines, bool parents = false)
{
    int i, delta;
    float pow, lambda = 1;
    if(kind_of_fine == 1)
        for (i = 0; i < a; i++)
            lambda*=C*globali;
    if(kind_of_fine == 2)
    {

    }
    if(kind_of_func == 0||kind_of_func == 1||kind_of_func == 2||kind_of_func == 3||kind_of_func == 4)
        delta = -1;
    if(kind_of_func == 5||kind_of_func == 6)
        delta = 1;
    if(parents)
    {
        for (i = 0; i < n; i++)
        {
            power[i] = f[i]+delta*lambda*fines[i];
        }
    }
    else
    {
        for (i = 0; i < n; i++)
        {
            power[i+n] = f[i+n]+delta*lambda*fines[i+n];
        }
    }

}

void limits(float*left, float *h, int kind_of_func, float* opt1, float* opt2, float e, int* number_m, int number_x, int *m=0)
{
    float nparts, z;
    int s = 1, i, temp_m;
    float *right = new float[number_x];
    if (kind_of_func == 0)
    {
        left[0] = 0;
        right[0] = 5;
        left[1] = -20;
        right[1] = 20;
        *opt1 = 3.65148371670110760;
        *opt2 = -6.6666666666666666;
        z = 69.8481705003444092;
    }
    if (kind_of_func == 1)
    {
        left[0] = 0;
        right[0] = 5;
        left[1] = 0;
        right[1] = 8;
        *opt1 = 4;
        *opt2 = 7.989358247;
        z = 79.82984520;
    }
    if (kind_of_func == 2)
    {
        left[0] = 0;
        right[0] = 3;
        left[1] = -1;
        right[1] = 6;
        *opt1 = 2.16666;
        *opt2 = 0.66666;
        z = 11.16666667;
    }
    if (kind_of_func == 3)
    {
        left[0] = 0;
        right[0] = 3;
        left[1] = -1;
        right[1] = 6;
        *opt1 = 2.16666;
        *opt2 = 0.66666;
        z = 11.16666667;
    }
    for (i = 0; i < number_x; i++)
    {
        nparts = (right[i] - left[i])/e;
        while (s < nparts + 1)
        {
            s *= 2;
            temp_m++;
        }
        number_m[i] = temp_m;
        h[i] = (right[i] - left[i])/(pow(2,temp_m)-1);
        s = 1;
        *m+=temp_m;
        temp_m = 0;
    }
    delete[] right;
}


int main()
{
    //ofstream fout("2.txt");

	setlocale(0, "");
	srand(time(NULL));
	string sel_switch = "tour";// prop, rang, tour
	int cross_switch = 3;//1, 2, 3
	string mut_switch = "strong";//weak, average, strong
	// 2 плохо
	int kind_of_func = 2;
	int kind_of_fine = 1; // 1-динамические, 2-адаптивные
	int n = 50, m = 0, i, j, in, k, fights = n*2, globali, i_launch, number_of_popul = 100, gener, counter, without_fines = 0;
	float e = 0.01,  h1, h2, opt1, opt2, number = 0, number_of_launch = 100, reability,min_diff1, min_diff2, x1min, x2min,
	powermin, sum_gener = 0, average_gener, C = 0.5, a = 2., b = 2., b1 = 1.1, b2 = 1.2, kk = 5;
	bool flag;
	int number_x;
	if (kind_of_func == 0||kind_of_func == 1||kind_of_func == 2||kind_of_func == 4||kind_of_func == 6)
        number_x = 2;
    if(kind_of_func == 3||kind_of_func == 5)
        number_x = 3;



    float *left = new float[number_x];
    float *h = new float[number_x];
    int *number_m = new int[number_x];
    limits(left, h, kind_of_func, &opt1, &opt2, e, number_m, number_x, &m);


    int** A1 = new int* [n];
    for (i = 0; i < n; i++)
        A1[i] = new int[m];
    int** A2 = new int* [n];
    for (i = 0; i < n; i++)
        A2[i] = new int[m];
    float *power1 = new float[n*2];
    float *fines = new float[n*2];
    float *fines2 = new float[n];
    //float *power2 = new float[n];
    int** parents = new int*[fights];
    for (i = 0; i < fights; i++)
        parents[i] = new int[m];
    int** children = new int* [fights];
    for (i = 0; i < fights; i++)
        children[i] = new int[m];
    int *index = new int[n*2];
    float** x = new float* [n*2];
    for (i = 0; i < n*2; i++)
        x[i] = new float[number_x];
    float** x2 = new float* [n];
    for (i = 0; i < n; i++)
        x2[i] = new float[number_x];
    float *f = new float[n*2];
    float *f2 = new float[n];




    //инициализация
    for (i = 0; i < n; i++)
        for (j = 0; j < m; j++)
            A1[i][j] = rand() % 2;

    count_x(A1, n, h, left, x, number_m, number_x, true);
    func(x, f, n, kind_of_func, number_x, true);
    fine(x, n, fines, b, kind_of_func, true);
    //цикл поколений
    for (globali = 1; globali < number_of_popul; globali++)
    {
        fitness(power1, f, globali, fines, C, a, b, n, kind_of_func, kind_of_fine, b1, b2, k, without_fines, true);
        cout << "Generation " << globali << endl;
        for(i = 0; i < n; i++)
        {
            for (j = 0; j < number_x; j++)
                cout << x[i][j] << "  ";
            cout << f[i] << "  " << fines[i] << "  " << power1[i] << endl;
        }
        selection (sel_switch, power1, fights, n, parents, A1, m);
        cross(fights, m, cross_switch, children, parents);
        mutation(mut_switch, n, m, children);
        count_x(children, n, h, left, x, number_m, number_x);
        func(x, f, n, kind_of_func, number_x);
        fine(x, n, fines, b, kind_of_func);
        fitness(power1, f, globali, fines, C, a, b, n, kind_of_func, kind_of_fine, b1, b2, k, without_fines);
        for (i = 0; i < n*2; i++)
            index[i] = i;
        quicksort(power1, index, 0, fights - 1);
        //создание новой популяции
        for (i = 0, k = fights-1; i < n; i++,k--)
        {
            in = index[k];
            if (in < n)
            {
                for (j = 0; j < m; j++)
                    A2[i][j] = A1[in][j];
            }
            if (in >= n)
            {
                 for (j = 0; j < m; j++)
                    A2[i][j] = children [in-n][j];
            }
            for (j = 0; j < number_x; j++)
                x2[i][j] = x[in][j];
            fines2[i] = fines[in];
            f2[i] = f[in];
        }
       for (i = 0; i < n; i++)
       {
           for (j = 0; j < number_x; j++)
                x[i][j] = x2[i][j];
           fines[i] = fines2[i];
           f[i] = f2[i];
           for (j=0; j < m; j++)
           {
               A1[i][j]=A2[i][j];
           }
       }
        if (fines[0] == 0)
            without_fines+=1;
        else
            without_fines = 0;

        /*min_diff1 = abs(x1[0]-opt1);
        min_diff2 = abs(x2[0]-opt2);
        for (i = 0; i < n; i++)
        {
            if(abs(x1[i]-opt1) <= e&&abs(x2[i]-opt2) <= e)
            {
                number+=1;
                flag = true;
                gener = globali;
                sum_gener += gener;
                break;
            }
            if (abs(x1[i]-opt1)<= min_diff1&&abs(x2[i]-opt2)<= min_diff2)
            {
                min_diff1 = abs(x1[i]-opt1);
                min_diff2 = abs(x2[i]-opt2);
                x1min = x1[i];
                x2min = x2[i];
                powermin = power1[i];

            }
        }*/
    }
    //reability = number/number_of_launch;
    //average_gener = ceil(sum_gener/number_of_launch);
    //cout << "Функция " <<kind_of_func<<" селекция "<<sel_switch<<" скрещивание "<<cross_switch<<" мутация "<<mut_switch<<endl;
    //cout << "Надежность "<<reability << endl;
    //fout << reability <<" "<< average_gener<< '\t';
    //cout << "В среднем решение находится на  " << average_gener << " поколении" << endl;

    //fout << endl;
    for (i = 0; i < n; i++)
        delete [] A2[i];
    delete[] A2;
    for (i = 0; i < fights; i++)
        delete[] children[i];
    delete[] children;
    for (i = 0; i < fights; i++)
        delete[] parents[i];
    delete[] parents;
    delete[] index;
    //delete[] power2;
    delete[] power1;
    for (i = 0; i < n; i++)
        delete [] A1[i];
    delete[] A1;
    delete[] fines;
    delete[] fines2;
    delete[] f;
    delete[] f2;
    delete[] left;
    delete[] h;
    delete[] number_m;
    for (i = 0; i < n; i++)
        delete [] x[i];
    delete[] x;
    for (i = 0; i < n; i++)
        delete [] x2[i];
    delete[] x2;

    //fout.close();
}
