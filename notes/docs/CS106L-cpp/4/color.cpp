// #define DEFINE_COLOR(color,opposite) case color: return opposite;
#include <iostream>
using namespace std;

enum Color {
  #define DEFINE_COLOR(color,opposite) color,
  #include "color.h"
  #undef DEFINE_COLOR
};
// enum Color {Red,Yellow,Green,Blue,Magenta,Cyan};
Color GetOppositeColor(Color c){
  switch (c){
    #define DEFINE_COLOR(color,opposite) case color: return opposite;
    #include "color.h"
    #undef DEFINE_COLOR
    default: return c;
  }
}
string ColorToString(Color c){
  switch(c) {
    #define DEFINE_COLOR(color,opposite) case color: return #color;
    #include "color.h";
    #undef DEFINE_COLOR;
    default: return "<unknown>";
  }
}
