#include <iostream>

#include "taggedarray.h"

using namespace std;

int main()
{
    TaggedArray<int, double, float> map {1, 2};

    const auto& cmap = map;

    cout << "double: " << cmap.get<double>() << "\n";
    cout << "float: " << map.get<float>() << "\n";
}
