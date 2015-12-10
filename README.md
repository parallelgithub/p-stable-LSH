# Introduction

本專案練習使用 p-stable distribution 實作 Locality-sensitive hashing，
主要是參考2004paper以及其名為 E2LSH 的實作，後者在實作時略為修改方法，
原paper是解approximation version R-NN problem，E2LSH則是解random version R-NN problem。
本專案主要依據的是E2LSH的實作方法。

本專案解m-nearest neighbors problem，因此E2LSH中的R參數...
E2LSH讓使用者輸入R值並訂R為目標，演算法用R與p算出機率lower bound，以求出欲達成目標所需的 k, L 值

程式碼根據預先給定的所有vector建立hash table，query vector再從table裡計算nearest neighbors。

建立hash table的過程大致分為兩步驟。第一步先將vector降低維度，第二步再將降維的vector計算hash value存進table。
第一步是降維用的方法是 Locality-sensitive hashing (LSH)，第二步就是一般的hashing。

p function是p-stable LSH實作的核心之一。
p function是由paper而來，在paper中有定義及其積分形式 在E2LSH中有實作
它提供了hash collision機率的lower bound，
因此也可幫我們推算k, L，
以保證整個LSH scheme成功的機率。


# Input Files format
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
*	Functional style vs. Imperative style

	When I implement the inner product of two vectors, I found that the functional style programming is slower than the imperative style programming.
functional coding :
		val dot = (vectorV zip x.vectorA).map{case (a,b) => a * b}.foldLeft(0.0)(_+_)
imperative coding :
		for(i <- 0 until len) sum = sum + v1(i) * v2(i)

	Functional coding style in scala sometimes may not be slow: http://stackoverflow.com/questions/2794823/is-scala-functional-programming-slower-than-traditional-coding

	"If we wanted to write the fastest code possible, we would essentially write Scala as if it were C."
	"The critical trade-off for us is that writing clean Scala is faster, less error prone, and easier to maintain than the alternative."
	https://www.sumologic.com/2012/07/23/3-tips-for-writing-performant-scala/

*	Read code of E2LSH
	https://www.facebook.com/notes/%E6%B4%AA%E5%A3%AB%E7%81%9D/%E6%82%85%E8%AE%80%E7%A8%8B%E5%BC%8F/1054461137918361