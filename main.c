/*
 新粒子群算法求解tsplib中的rand400
 屏蔽掉了原粒子群算法中的速度，改用随机3种不同的操作进行迭代进化
 1.拷贝pbest为自身，进行变异操作
 2.拷贝gbest为自身，进行变异操作
 3.直接进行变异操作
 配合邻域2-opt优化，最终结果与已知最优解误差在0.5%左右
 这种新型粒子群算法大大降低了编程难度，尤其对解决复杂离散问题效果非常好
 如有疑问请联系 QQ281513754
 */
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define TEST_NUM    10          //实验次数
#define USE_2OPT    0           //是否使用2opt优化 1使用0不使用
#define PRINT_PER_LOOP (USE_2OPT ? 10 : 100) //每隔多少代打印一次结果
#define DATA_PATH   "/Users/zhixiang/Desktop/PSO_TSP_C/rand400.txt"  //读取城市坐标文件的路径
#define RANDOM01    (rand() % 32768 / 32768.0) //[0,1)间的随机小数
#define CITY_NUM    400         //城市数
#define P_NUM   1000            //粒子个数
#define C1  0.3                 //参数
#define C2  0.6                 //参数
#define SAME_STOP   (USE_2OPT ? 200 : 1000)        // 连续多少代最优解相同时终止计算
#define ERROR_STOP  0.000001    //误差率小于多少时终止计算
#define MAX_ITERATION  10000     //最多迭代次数
#define RANK_MAX    20          //2opt邻域搜索点个数
#define BEST_DISTANT    14722.0
int _cityPoint[CITY_NUM][2];
int _gbestIndex;                        //gbest下标
double _fitness[P_NUM];                 //每个粒子pbest的目标函数值，即路程
int _pos[P_NUM][CITY_NUM];              //当前粒子城市序列
int _pest[P_NUM][CITY_NUM];             //粒子的pbest序列
double _cityDis[CITY_NUM][CITY_NUM];    //城市间距离矩阵
int _rank[CITY_NUM][RANK_MAX];          //离某个城市最近的若干个城市排序
//读取城市坐标数据
void readCityData() {
    FILE *fp;
    if ((fp = fopen(DATA_PATH, "r")) != NULL)
    {
        char buf[64] = "";
        int count = 0;
        while (fscanf(fp, "%s", buf) != EOF) {
            if (count % 3 == 1) {
                _cityPoint[count / 3][0] = atoi(buf);
            } else if (count % 3 == 2) {
                _cityPoint[count / 3][1] = atoi(buf);
            }
            count++;
        }
    }
    else {
        printf("open txt failed\n");
        exit(-1);
    }
    fclose(fp);
}
//计算目标函数
double calFitness(int tsp[CITY_NUM]){
    double total = 0;
    for (int i = 0; i < CITY_NUM - 1; i++) {
        total += _cityDis[tsp[i]][tsp[i + 1]];
    }
    total += _cityDis[tsp[0]][tsp[CITY_NUM - 1]];
    return total;
}
//变异
void mutate(int p[CITY_NUM]) {
    int r = rand() % 6;
    int index1, index2;
    do{
        index1 = rand() % CITY_NUM;
        index2 = rand() % CITY_NUM;
    } while (index1 >= index2);
    int temp;
    if (r == 0) {   //插入变异1
        temp = p[index1];
        memmove(p + index1, p + index1 + 1, sizeof(int) * (index2 - index1));
        p[index2] = temp;
    } else if (r == 1) {   //插入变异2
        temp = p[index2];
        memmove(p + index1 + 1, p + index1, sizeof(int) * (index2 - index1));
        p[index1] = temp;
    } else if (r == 2) {   //交换变异
        temp = p[index1];
        p[index1] = p[index2];
        p[index2] = temp;
    } else if (r == 3) {   //逆序变异
        for (int i = index1; i <= (index1 + index2) / 2; i++) {
            temp = p[i];
            p[i] = p[index2 + index1 - i];
            p[index2 + index1 - i] = temp;
        }
    } else if (r == 4) {    //连续乱序变异
        int len =  rand() % (CITY_NUM / 10);
        len = 2 + rand() % (len + 1);
        index2 = index1 + len;
        if (index2 >= CITY_NUM) {
            index2 = CITY_NUM - 1;
        }
        for (int i = index1; i < index2; i++) {
            int pos = i + rand() % (index2 + 1 - i);
            temp = p[i];
            p[i] = p[pos];
            p[pos] = temp;
        }
    }  else if (r == 5) { //区域乱序变异
        int series[CITY_NUM] = { 0 };
        int seriesNum = 0;
        double r = 20 + rand() % 300;
        series[seriesNum++] = index1;
        for (int i = 0; i < CITY_NUM; i++) {
            if (_cityDis[i][index1] < r) {
                series[seriesNum++] = i;
            }
        }
        for (int i = 0; i < seriesNum - 1; i++) {
            int pos = i + rand() % (seriesNum - i);
            temp = p[series[i]];
            p[series[i]] = p[series[pos]];
            p[series[pos]] = temp;
        }
    }
}
//更新适应值
void updateFitness(){
    for (int i = 0; i < P_NUM; i++){
        double fit = calFitness(_pos[i]);
        if (fit < _fitness[i] || (_gbestIndex != i && RANDOM01 < 0.01)){
            _fitness[i] = fit;
            memcpy(_pest[i], _pos[i], sizeof(int) * CITY_NUM);
        }
        if (fit < _fitness[_gbestIndex]){
            _gbestIndex = i;
        }
    }
}
//初始化粒子
void initRes() {
    _gbestIndex = 0;
    for (int i = 0; i < P_NUM; i++) {
        _fitness[i] = 1e12;
    }
    for (int n = 0; n < P_NUM; n++) {
        for (int i = 0; i < CITY_NUM; i++) {
            _pos[n][i] = i;
        }
        for (int i = CITY_NUM - 1; i > 0; i--) {
            int index = rand() % (i + 1);
            int temp = _pos[n][index];
            _pos[n][index] = _pos[n][i];
            _pos[n][i] = temp;
        }
        memcpy(_pest[n], _pos[n], sizeof(int) * CITY_NUM);
    }
}
//2_opt优化
void search_2opt(int *tsp) {
    int check[CITY_NUM] = {0};
    int indexs[CITY_NUM];
    for (int k = 0; k <CITY_NUM; k++) {
        indexs[tsp[k]] = k;
    }
    int flag;
    int count = 0;
    do {
        count++;
        flag = 0;
        for (int p = 0; p <CITY_NUM; p++) {
            int city_i = tsp[p];
            if (check[city_i]) {
                continue;
            }
            int pos_i = p;
            int city_aft_i = tsp[(p + 1) % CITY_NUM];
            double search_range_i = _cityDis[city_i][city_aft_i];
            for (int q = 0; q < RANK_MAX; q++) {
                int city_j = _rank[city_i][q];
                if (check[city_j]) {
                    continue;
                }
                if (_cityDis[city_i][city_j] > search_range_i - 0.0000001 || city_j == city_aft_i) {
                    break;
                }
                int pos_j = indexs[city_j];
                int city_aft_j = tsp[(pos_j + 1) % CITY_NUM];
                double before = _cityDis[city_i][city_aft_i] + _cityDis[city_j][city_aft_j];
                double after = _cityDis[city_i][city_j] + _cityDis[city_aft_i][city_aft_j];
                if (after < before) {
                    int i = pos_i;
                    int j = pos_j;
                    flag = 1;
                    if (i > j) {
                        int swap = j;
                        j = i;
                        i = swap;
                    }
                    check[city_i] = 0;
                    check[city_aft_i] = 0;
                    check[city_j] = 0;
                    check[city_aft_j] = 0;
                    int m = i + 1;
                    int n = j;
                    while (m < n) {
                        int t1;
                        int t2;
                        t1 = tsp[m];
                        t2 = tsp[n];
                        tsp[m] = t2;
                        tsp[n] = t1;
                        indexs[t1] = n;
                        indexs[t2] = m;
                        m++;
                        n--;
                    }
                    city_aft_i = city_j;
                    search_range_i = _cityDis[city_i][city_j];
                }
            }
            if (!flag) {
                check[city_i] = 1;
            }
        }
    } while (flag);
}
//每次迭代之后对每一个解进行相关概率的改变
void NPSO() {
    for (int i = 0; i < P_NUM; i++) {
        double r = RANDOM01;
        if (r < C1) {
            memcpy(_pos[i], _pest[i], sizeof(int) * CITY_NUM);
        } else if (r < C2){
            memcpy(_pos[i], _pest[_gbestIndex], sizeof(int) * CITY_NUM);
        }
        mutate(_pos[i]);
        if (USE_2OPT) {
            search_2opt(_pos[i]);
        }
    }
}
//初始化距离及邻接矩阵
void initDis() {
    memset(_cityDis, 0, sizeof(double) * CITY_NUM * CITY_NUM);
    for (int i = 0; i <CITY_NUM; i++) {
        for (int j = 0; j < i; j++) {
            double dx = _cityPoint[i][0] - _cityPoint[j][0];
            double dy = _cityPoint[i][1] - _cityPoint[j][1];
            _cityDis[i][j] = _cityDis[j][i] = sqrt(dx * dx + dy * dy);
        }
        _cityDis[i][i] = 1e10;
    }
    memset(_rank, 0, sizeof(int) * CITY_NUM * RANK_MAX);
    int series[CITY_NUM];
    double tempDis[CITY_NUM] = { 0 };
    for (int i = 0; i < CITY_NUM; i++) {
        for (int j = 0; j < CITY_NUM; j++) {
            series[j] = j;
        }
        memcpy(tempDis, _cityDis[i], sizeof(double) *CITY_NUM);
        for (int j = 0; j < RANK_MAX; j++) {
            for (int k = j + 1; k <CITY_NUM; k++) {
                if (tempDis[k] < tempDis[j]) {
                    int temp = series[j];
                    series[j] = series[k];
                    series[k] = temp;
                    double tempD = tempDis[j];
                    tempDis[j] = tempDis[k];
                    tempDis[k] = tempD;
                }
            }
        }
        memcpy(_rank[i], series, sizeof(int) * RANK_MAX);
    }
}
int main(int argc, const char * argv[]){
    srand((unsigned int)time(NULL));
    readCityData();
    initDis();
    int loopNum[TEST_NUM];
    double bestFit[TEST_NUM];
    double error[TEST_NUM];
    double timeUsed[TEST_NUM];
    for (int t = 0; t < TEST_NUM; t++) {
        initRes();
        struct timeval startTime;
        gettimeofday(&startTime, NULL);
        double lastBest = 1e10;
        int same = 0;
        for (int n = 0; n < MAX_ITERATION; n++){
            NPSO();
            updateFitness();
            if (lastBest - _fitness[_gbestIndex] < 1e-10) {
                same++;
            } else {
                same = 0;
            }
            lastBest = _fitness[_gbestIndex];
            int end = 0;
            if (same >= SAME_STOP || (_fitness[_gbestIndex] - BEST_DISTANT) * 100 / BEST_DISTANT < ERROR_STOP){
                end = 1;
            }
            if (n % PRINT_PER_LOOP == PRINT_PER_LOOP - 1 || end) {
                struct timeval endTime;
                gettimeofday(&endTime, NULL);
                double runningTime = endTime.tv_sec - startTime.tv_sec + (endTime.tv_usec - startTime.tv_usec) / 1000000.0;
                printf("%d-%4d %.2f, %.3lf\n", t + 1, n + 1, _fitness[_gbestIndex], runningTime);
                timeUsed[t] = runningTime;
                loopNum[t] = n + 1;
            }
            if (end) {
                break;
            }
        }
        bestFit[t] = _fitness[_gbestIndex];
        error[t] = (_fitness[_gbestIndex] - BEST_DISTANT) / BEST_DISTANT;
    }
    printf("实验\t迭代次数\t路径长度\t与最优解误差\t用时\n");
    for (int t = 0; t < TEST_NUM; t++) {
        printf("%d\t%d\t%.2f\t%.2f%%\t%.2f秒\n", t + 1, loopNum[t], bestFit[t], error[t] * 100, timeUsed[t]);
    }
    double avgLoopNum = 0;
    double avgBestFit = 0;
    double avgError = 0;
    double avgTimeUsed = 0;
    for (int t = 0; t < TEST_NUM; t++) {
        avgLoopNum += loopNum[t] * 1.0/ TEST_NUM;
        avgBestFit += bestFit[t] / TEST_NUM;
        avgError += error[t] / TEST_NUM;
        avgTimeUsed += timeUsed[t] / TEST_NUM;
    }
    printf("平均\t%.f\t%.2f\t%.2f%%\t%.2f秒\n",  avgLoopNum, avgBestFit, avgError * 100, avgTimeUsed);
    return 0;
}
