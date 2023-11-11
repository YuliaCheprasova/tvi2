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



void count_x(int **A, int n, int m, float h1, float h2, float left1, float left2, float *x1, float *x2, int m1, int m2, bool parents = false)
{
    int temp = 0, i, j, k;
    for (i = 0; i < n; i++)
    {
        for (j = 0, k = m1-1; j < m1; j++, k--)
        {
            if (A[i][j] == 1)
            {
                temp += pow(2,k);
            }
        }
        if (parents)
            x1[i] = h1*temp + left1;
        else
            x1[i+n] = h1*temp + left1;
        temp = 0;
        for (j = m1, k = m2-1; j < m; j++, k--)
        {
            if (A[i][j] == 1)
            {
                temp += pow(2,k);
            }
        }
        if (parents)
            x2[i] = h2*temp + left2;
        else
            x2[i+n] = h2*temp + left2;
        temp = 0;
    }
}

void func(float *x1, float *x2, float *f, int n, int kind_of_func, bool parents = false)
{
    int i, j = 0;
    if (parents == false)
    {
        j = n;
        n = n*2;
    }
    if (kind_of_func == 0)
        for (i = j; i < n; i++)
            f[i] = 10*x1[i]-5*x2[i];
}

void fine(float *x1, float *x2, int n, float *fines, float b, int kind_of_func, bool parents = false)
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
                mult_fines*=max(float(0), x2[i]-15);
            temp_fines+=mult_fines;
            mult_fines = 1;
            for (j = 0; j < b; j++)
                mult_fines*=max(float(0), x2[i]+2*x1[i]*x1[i]-20);
            temp_fines+=mult_fines;
            mult_fines = 1;
            for (j = 0; j < b; j++)
                mult_fines*=max(0., -((x1[i]*x1[i])/2.)-x2[i]);
            temp_fines+=mult_fines;
            fines[i] = temp_fines;
        }

    }
}

void fitness(float *power, float* f, int globali, float *fines, float C, float a, float b, int n, int kind_of_func, bool parents = false)
{
    int i, delta;
    float pow, lambda = 1;
    for (i = 0; i < a; i++)
    {
        lambda*=C*globali;
    }
    if(kind_of_func == 0)
        delta = -1;
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

void limits(float*left1, float *left2, float *h1, float *h2, int kind_of_func, float* opt1, float* opt2, float e, int* m1=0, int* m2=0, int *m=0)
{
    float nparts, right1, right2, z;
    int s = 1;
    if (kind_of_func == 0||kind_of_func == 1||kind_of_func == 3)
    {
        *left1 = 0;
        right1 = 5;
        *left2 = -20;
        right2 = 20;
        *opt1 = 3.65148371670110760;
        *opt2 = -6.6666666666666666;
        z = 69.8481705003444092;
    }
    /*if (kind_of_func == 2)
    {
        *left = -2;
        *right = 2;
        *opt1 = 1;
        *opt2 = 1;
    }
    if (kind_of_func == 4)
    {
        *left = -10;
        *right = 10;
        *opt1 = 1;
        *opt2 = 1;
    }
    if (kind_of_func == 5)
    {
        *left = -10;
        *right = 10;
        *opt1 = 0;
        *opt2 = 0;
    }
    if (kind_of_func == 6)
    {
        *left = -2.5;
        *right = 2.5;
        *opt1 = 0;
        *opt2 = 0;
    }
    if (kind_of_func == 7)
    {
        *left = -5;
        *right = 5;
        *opt1 = 0;
        *opt2 = 0;
    }
    if (kind_of_func == 8)
    {
        *left = -4;
        *right = 4;
        *opt1 = 0;
        *opt2 = 0;
    }
    if (kind_of_func == 9)
    {
        *left = 0;
        *right = 4;
        *opt1 = 0;
        *opt2 = 0;
    }
    if (kind_of_func == 10||kind_of_func == 11)
    {
        *left = 0;
        *right = 4;
        *opt1 = 1.99516;
        *opt2 = 1.99516;
    }
    if (kind_of_func == 12)
    {
        *left = -65;
        *right = 65;
        *opt1 = -32;
        *opt2 = -32;
    }*/

    nparts = (right1 - *left1)/e;
    while (s < nparts + 1)
    {
        s *= 2;
        (*m1)++;
    }
    *h1 = (right1 - *left1)/(pow(2,*m1)-1);
    s = 1;
    nparts = (right2 - *left2)/e;
    while (s < nparts + 1)
    {
        s *= 2;
        (*m2)++;
    }
    *h2 = (right2 - *left2)/(pow(2,*m2)-1);
    *m = *m1 + *m2;
}


int main()
{
    //ofstream fout("2.txt");

	setlocale(0, "");
	srand(time(NULL));
	string sel_switch = "tour";// prop, rang, tour
	int cross_switch = 3;//1, 2, 3
	string mut_switch = "strong";//weak, average, strong
	int kind_of_func = 0;
	int n = 50, m = 0, m1 = 0, m2 = 0, i, j, in, k, fights = n*2, globali, i_launch, number_of_popul = 100, gener, counter;
	float left1, left2, e = 0.001,  h1, h2, opt1, opt2, number = 0, number_of_launch = 100, reability,
	min_diff1, min_diff2, x1min, x2min, powermin, sum_gener = 0, average_gener, C = 0.5, a = 2., b = 2.;
	bool flag;



    limits(&left1, &left2, &h1, &h2, kind_of_func, &opt1, &opt2, e, &m1, &m2, &m);


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
    int *index = new int[n + fights];
    float *x1 = new float[n*2];
    float *x2 = new float[n*2];
    float *x12 = new float[n];
    float *x22 = new float[n];
    float *f = new float[n*2];
    float *f2 = new float[n];




    //инициализаци€
    for (i = 0; i < n; i++)
        for (j = 0; j < m; j++)
            A1[i][j] = rand() % 2;

    count_x(A1, n, m, h1, h2, left1, left2, x1, x2, m1, m2, true);
    func(x1, x2, f, n, kind_of_func, true);
    fine(x1, x2, n, fines, b, kind_of_func, true);
    //цикл поколений
    for (globali = 1; globali < number_of_popul; globali++)
    {
        fitness(power1, f, globali, fines, C, a, b, n, kind_of_func, true);
        cout << "Generation " << globali << endl;
        for(i = 0; i < n; i++)
        {
            cout << x1[i] << " " << x2[i] << " " << f[i] << " " << fines[i] << " " << power1[i] << endl;
        }
        selection (sel_switch, power1, fights, n, parents, A1, m);
        cross(fights, m, cross_switch, children, parents);
        mutation(mut_switch, n, m, children);
        count_x(children, n, m, h1, h2, left1, left2, x1, x2, m1, m2);
        func(x1, x2, f, n, kind_of_func);
        fine(x1, x2, n, fines, b, kind_of_func);
        fitness(power1, f, globali, fines, C, a, b, n, kind_of_func);
        for (i = 0; i < n*2; i++)
            index[i] = i;
        quicksort(power1, index, 0, fights - 1);
        //создание новой попул€ции
        for (i = 0, k = fights-1; i < n; i++,k--)
        {
            in = index[k];
            if (in < n)
            {
                for (j = 0; j < m; j++)
                    A2[i][j] = A1[in][j];
                x12[i] = x1[in];
                x22[i] = x2[in];
                fines2[i] = fines[in];
                f2[i] = f[in];
            }
            if (in >= n)
            {
                 for (j = 0; j < m; j++)
                    A2[i][j] = children [in-n][j];
                x12[i] = x1[in];
                x22[i] = x2[in];
                fines2[i] = fines[in];
                f2[i] = f[in];
            }
        }
       for (i = 0; i < n; i++)
       {
           x1[i] = x12[i];
           x2[i] = x22[i];
           fines[i] = fines2[i];
           f[i] = f2[i];
           for (j=0; j < m; j++)
           {
               A1[i][j]=A2[i][j];
           }
       }


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
    //cout << "‘ункци€ " <<kind_of_func<<" селекци€ "<<sel_switch<<" скрещивание "<<cross_switch<<" мутаци€ "<<mut_switch<<endl;
    //cout << "Ќадежность "<<reability << endl;
    //fout << reability <<" "<< average_gener<< '\t';
    //cout << "¬ среднем решение находитс€ на  " << average_gener << " поколении" << endl;

    //fout << endl;
    delete[] x1;
    delete[] x2;
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
    delete[] x12;
    delete[] x22;
    delete[] fines;
    delete[] fines2;
    delete[] f;
    delete[] f2;

    //fout.close();
}
// разобратьс€ почему получаютс€ нули
