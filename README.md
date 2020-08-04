# SimpleNeuralNetwork
This code is developed based on the code by huangzehou: https://github.com/huangzehao/SimpleNeuralNetwork which is the implementation of simple neural network based on video [Neural Net in C++ Tutorial](https://vimeo.com/19569529) by David Miller. 
# Test in Linux
1 Gernerate training data to slove XOR problem
```
    root -l write_text.C
```
2 Test neural netwrok
```
    g++ ./adam-test.cpp -o testann2 `root-config --cflags --glibs`
    ./testann2
```
And you will get the result!

