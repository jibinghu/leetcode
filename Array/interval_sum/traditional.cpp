#include <iostream>
#include <vector>
// O(n*m)
int main()
{
    int n, m;
    std::cin >> n;
    std::vector<int> array;
    for (int i = 0; i < n; i++)
    {
        std::cin >> m;
        array.push_back(m);
    }
    int begin, end;
    while (std::cin >> begin >> end)
    {
        int sum = 0;
        for (int i = begin; i < end + 1; i++)
        {
            sum += array.at(i);
        }
        std::cout << sum << std::endl;
    }
}