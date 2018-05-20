clear all
mata: 

// void test(in) 
// {
// 
// if (something) class myclass1 scalar C
// else           class myclass2 scalar C
// 
// if (something) real scalar x
// else           real matrix x
//  
//  x = 2
//  x
// }

class myclass1 {
  void new()
  real scalar x
}
class myclass2 extends myclass1 {
  void new()
  real scalar x
}
void myclass1::new() {
  "hello from 1"
}
void myclass2::new() {
  "hello from 2"
}

void test1() {
  "start"
  class myclass1 scalar C1
  class myclass2 scalar C2
  pointer C
  "declarations done"
  if (0) C = &C1
  else   C = &C2
//  C2.x = 2
////  *C.x = 2
//  C->x = 2

  "end"
}
"run test 1"
test1() 

void test2() {
  if (0) class myclass1 scalar C
  else   class myclass2 scalar C
}
"run test 2"
test2()

end

