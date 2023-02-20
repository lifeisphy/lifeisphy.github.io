#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

const int NUM_LINES = 4; // how to read arbitary lines?
const int NUM_COLUMNS = 3;
const int COLUMN_WIDTH = 20;
void PrintTableBody()
{
  ifstream input("data");
  // for(int k=0; k<NUM_LINES; ++k){
  //   int INTVAL;
  //   double DOUBLEVAL;
  //   input>>INTVAL>>DOUBLEVAL;
  //   cout<<setw(COLUMN_WIDTH)<<(k+1)<<" | ";
  //   cout<<setw(COLUMN_WIDTH)<<INTVAL<<" | ";
  //   cout<<setw(COLUMN_WIDTH)<<DOUBLEVAL<<endl;
  // }
  // arbitary lines:
  // int lineno= 1;
  // while (true)
  // {
  //   int intval;
  //   double doubleval;
  //   input >> intval >> doubleval;
  //   if (input.fail())
  //     break;
  //   cout << setw(COLUMN_WIDTH) << lineno++ << " | ";
  //   cout << setw(COLUMN_WIDTH) << intval << " | ";
  //   cout << setw(COLUMN_WIDTH) << doubleval << endl;
  // }
  // or:
  int intval;
  double doubleval;
  while (input>>intval>>doubleval) // return false when encounters failure
  {
    if (input.fail())
      break;
    cout << setw(COLUMN_WIDTH) << lineno++ << " | ";
    cout << setw(COLUMN_WIDTH) << intval << " | ";
    cout << setw(COLUMN_WIDTH) << doubleval << endl;
  }
}
void PrintTableHeader()
{
  // for(int col=0;col<NUM_COLUMNS;++col){
  //   for(int k=0; k<COLUMN_WIDTH; ++k){
  //     cout<< '-';
  //   }
  //   if(col!= NUM_COLUMNS-1){
  //     cout<<"-+-";
  //   }
  // }
  // or:
  for (int col = 0; col < NUM_COLUMNS; col++)
  {
    cout << setfill('-') << setw(COLUMN_WIDTH);
    cout << "";
    if (col != NUM_COLUMNS - 1)
      cout << "-+-";
  }
  cout << setfill(' ');
  cout << endl;
}
int main()
{
  PrintTableHeader();
  PrintTableBody();
}