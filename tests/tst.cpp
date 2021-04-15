#include <iostream>

#include "taggedtuple.h"

using namespace std;

int main()
{
    TaggedTuple<int, double, float> map(1, 2);

    const auto& cmap = map;

    cout << "double: " << get<double>(cmap) << ", float: " << get<float>(map) << "\n";
}
