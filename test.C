

int test1(){return 1;}

template <typename T>
void printVector(const std::vector<T>& vector_){
    std::cout << "{ ";
    bool isFirst = true;
    for(const auto& element: vector_){
      if(not isFirst) std::cout << ", ";
      else isFirst = false;
      std::cout << element;
    }
    std::cout << " }" << std::endl;
  }

void test(){

}
