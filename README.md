# Introduction

本專案練習使用 p-stable distribution 實作 Locality-sensitive hashing，簡稱LSH。

LSH演算法用來解決 nearest neighbor problem，這個 problem 在許多領域皆有廣泛的應用：如 data compression, databases & data mining, information retrieval, image & video databases, machine learning, pattern recognition, statistics and data analysis。LSH演算法的價值在於當 nearest neighbor problem 為高維度時，提供了 sublinear-time 的解法，因此LSH適合用於高維度的特徵處理。

本專案主要依據 E2LSH 的實作方法: 
http://www.mit.edu/~andoni/LSH/
我們參考了 E2LSH 的原始碼及其 manual。

解 nearest neighbor problem 的流程大致如下：
根據給定的所有 vector 建立 hash table，然後每個 query vector 再從 hash table 裡尋找 nearest neighbors。LSH演算法是整個流程的關鍵步驟。

建立 hash table 的過程分為兩步驟：第一步先將vector降低維度，第二步再將降維的vector計算hash value存進table。第一步是降維用的方法就是 Locality-sensitive hashing (LSH)，第二步則是一般的hashing。

# Input Files format
There are two input files. One contains base vectors, and the other contains query vectors.
This project find the nearest neighbor of each query vectors among all base vectors.

Both files have the same format.
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
1.	Functional style vs. Imperative style

	When I implement the inner product of two vectors, I found that the functional style programming is slower than the imperative style programming.
	* functional coding :
        `val dot = (vectorV zip x.vectorA).map{case (a,b) => a * b}.foldLeft(0.0)(_+_)`
	* imperative coding :
		`for(i <- 0 until len) sum = sum + v1(i) * v2(i)`

	Functional coding style in scala sometimes may not be slow: http://stackoverflow.com/questions/2794823/is-scala-functional-programming-slower-than-traditional-coding

	"If we wanted to write the fastest code possible, we would essentially write Scala as if it were C."
	"The critical trade-off for us is that writing clean Scala is faster, less error prone, and easier to maintain than the alternative."
	https://www.sumologic.com/2012/07/23/3-tips-for-writing-performant-scala/

	We compare several different ways of implementation of vector inner product: 
	https://github.com/parallelgithub/vector-inner-product

2.	Read code of E2LSH
	https://www.facebook.com/notes/%E6%B4%AA%E5%A3%AB%E7%81%9D/%E6%82%85%E8%AE%80%E7%A8%8B%E5%BC%8F/1054461137918361