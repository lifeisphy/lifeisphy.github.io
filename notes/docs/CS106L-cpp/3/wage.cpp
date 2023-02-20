#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <map>
#include <fstream>

using namespace std;
int main(int argc,char** argv){ 
  // try to input 2.718
  int age;
  double hourlyWage;
  cout<< "Please enter your wage: ";
  cin>>age;
  cout<< "Please enter your hourly wage: ";
  cin>>hourlyWage;
  cout<<age<<","<<hourlyWage;
  return 0;
}