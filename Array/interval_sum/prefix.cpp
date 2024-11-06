#include <vector>
#include <iostream>
// scanf 和 printf 在处理大量数据时通常优于 cin 和 cout
// O(n)
int main()
{
    int n, sum = 0;
    // std::cin >> n;
    scanf("%d", &n);
    std::vector<int> array(n), prefix(n);
    for (int i = 0; i < n; i++)
    {
        // std::cin >> array.at(i);
        scanf("%d", &array.at(i));
        sum += array.at(i);
        prefix.at(i) = sum;
    }
    int begin, end;
    while (std::cin >> begin >> end)
    {
        if (begin == 0)
        {
            // std::cout << prefix[end] << std::endl;
            printf("%d", prefix[end]);
        }
        else
        {
            // std::cout << prefix[end] - prefix[begin-1] << std::endl;
            printf("%d", prefix[end] - prefix[begin - 1]);
        }
        printf("\n");
    }
}