#include <iostream>
#include "dynamic/dynamic.hpp"

using namespace dyn;

int main(){
    std::cout << "Hello world!\n" ;
    packed_spsi my_spsi;
    my_spsi.push_back(10);
    my_spsi.push_back(20);
    my_spsi.push_back(30);

    std::cout << my_spsi.at(2) << "\n";
    return (EXIT_SUCCESS);
}
