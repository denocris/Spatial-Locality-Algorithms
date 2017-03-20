#include <string>
#include <iostream>
#include <bitset>

typedef std::bitset<64> bittype; // Alias

bool mybitf(const bittype & n, const unsigned int pos){
  return (n[pos] == 1);
}

void mybitset(bittype & n, const unsigned int pos){
  n[pos] = 1;
}

void mybitcrl(bittype & n, const unsigned int pos){
  n[pos] = 0;
}

void mybits(bittype & n, const unsigned int pos, const unsigned int lbit){
  bittype ret;
  for (int i = 0; i <= lbit; i ++){
    ret[i] = n[i + pos];
  }
  n = ret;
}

void mybitshiftR(bittype & n, const unsigned int len){
  bittype ret;
  for (int i = 0; i <= len + 1; i ++){
    ret[i + len] = n[i];
  }
  n = ret;
}

void mybitshiftL(bittype & n, const unsigned int len){
  bittype ret;
  for (int i = 0; i <= len + 1; i ++){
    ret[i] = n[i + len];
  }
  n = ret;
}



int main(){

  bittype mybit; // by defoult 0

  std::cout << mybit << std::endl;

  mybitset(mybit, 2);
  mybitset(mybit, 1);
  mybitset(mybit, 0);


  std::cout << mybit << std::endl;
  std::cout << mybit.to_ulong() << std::endl;

  mybitshiftR(mybit, 2);
  std::cout << mybit << std::endl;
  std::cout << mybit.to_ulong() << std::endl;

  // mybitshiftL(mybit, 2);
  // std::cout << mybit << std::endl;
  // std::cout << mybit.to_ulong() << std::endl;

  mybitcrl(mybit, 3);
  std::cout << mybit << std::endl;
  std::cout << mybit.to_ulong() << std::endl;

  mybits(mybit, 2, 3);
  std::cout << mybit << std::endl;
  std::cout << mybit.to_ulong() << std::endl;

return 0;

}
