#include <iostream>
#include <sstream>
using namespace std;

string GetAboutInformation(){
  stringstream result;
  result<<"This program was compiled on "<<__DATE__;
  result<<" at time "<<__TIME__;
  return result.str();
}
string GetAboutInformation2(){
  stringstream result;
  result<<__LINE__<<" "<<__FILE__;
  return result.str();
}
int main(){
  cout<<GetAboutInformation()<<endl;
  cout<<GetAboutInformation2()<<endl;
  return 0;
}