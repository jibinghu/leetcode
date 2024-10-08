#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>

int main() {
    /*
    蒙特卡洛算法估算圆周率π。假设在一个正方形内画一个内切圆，我们随机生成点，看这些点落在圆内的概率与正方形内所有点的比例关系，就可以估算出π的值
    */
    srand(time(0));  // 初始化随机数种子
    int total_points = 1000000;  // 总点数
    int points_in_circle = 0;    // 落在圆内的点数

    for (int i = 0; i < total_points; ++i) {
        double x = (double)rand() / RAND_MAX;  // 生成0到1之间的随机数
        double y = (double)rand() / RAND_MAX;

        if (x * x + y * y <= 1) {  // 判断是否在圆内
            points_in_circle++;
        }
    }

    double pi_estimate = 4.0 * points_in_circle / total_points;  // π的估计值
    std::cout << "Estimated value of pi: " << pi_estimate << std::endl;

    return 0;
}