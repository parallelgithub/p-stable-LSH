# Introduction

本專案練習使用 p-stable distribution 實作 Locality-sensitive hashing，
主要是參考2004paper以及其名為 E2LSH 的實作，後者在實作時略為修改方法，
原paper是approximation algorithm，E2LSH則實作為random algorithm。
本專案主要依據的是E2LSH的實作方法。

程式碼根據預先給定的所有vector建立hash table，query vector再從table裡計算nearest neighbors。

建立hash table的過程大致分為兩步驟。第一步先將vector降低維度，第二步再將降維的vector計算hash value存進table。
第一步是降維用的方法是 Locality-sensitive hashing (LSH)，第二步就是一般的hashing。

p function是p-stable LSH實作的核心之一。
p function是由paper而來，在paper中有定義及其積分形式 在E2LSH中有實作
它提供了hash collision機率的lower bound，
因此也可幫我們推算k, L，
以保證整個LSH scheme成功的機率。


# Files format
## Data set file
The file contains N lines, where N is the number of vertors. 
Each line represents the coordinate of a vertor.
A line contains M real numbers separated by a space, 
where M is the dimension of a vector.
The file format is as the following:
coordinate_1_of_point_1 coordinate_2_of_point_1 ... coordinate_D_of_point_1
coordinate_1_of_point_2 coordinate_2_of_point_2 ... coordinate_D_of_point_2
...
coordinate_1_of_point_N coordinate_2_of_point_N ... coordinate_D_of_point_N

# What I learn
1. functional coding style in scala sometimes may not be a good choise:
http://stackoverflow.com/questions/2794823/is-scala-functional-programming-slower-than-traditional-coding

functional coding :
val dot = (vectorV zip x.vectorA).map{case (a,b) => a * b}.foldLeft(0.0)(_+_)
imperative coding :
for(i <- 0 until len) sum = sum + v1(i) * v2(i)
	

number of inner-product operations is L*k and vector size is d, where L = 513, k = 45, d = 128

time consuming of two coding style:
188.2s / 16.2s = 11.6 times

2. read code of E2LSH