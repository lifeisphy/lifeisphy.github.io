#include <iostream>
using namespace std;
void printBinary(int a){
  for(int t=sizeof(int)*8-1;t>=0;t--){
    cout<<((a>>t)&1);
  }
  cout<<endl;
}
int bitManipulation2(int n, int i) {
  int a=~ n & (-1<<i);
  // int b=(!((n>>i)&1)<<i);
  int c=(n&~(-1<<i));
  printBinary(a);printBinary(c);
  printBinary(a|c);
  return (~ n & (-1<<(32-i)))|(n&~(-1<<(32-i)));
}

int main(int argc, char **argv) {
  printBinary(atoi(argv[1]));
  cout<< bitManipulation2(atoi(argv[1]),32-atoi(argv[2]))<<endl;
	return 0;
}