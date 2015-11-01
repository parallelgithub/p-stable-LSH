import scala.io.Source
import scala.util.Random
import scala.math


object LSH {

	val prime = 2147483647

	val binWidthW = 4.0
	val radiusR = 1.0
	val successProbability = 0.9

	//from http://picomath.org/scala/Erf.scala.html
	//   & http://www.johndcook.com/blog/cpp_erf/
	//It is correct, can be better...
	def ERF(x: Double): Double =  {
		val a1: Double =  0.254829592;
		val a2: Double = -0.284496736;
		val a3: Double =  1.421413741;
		val a4: Double = -1.453152027;
		val a5: Double =  1.061405429;
		val p: Double  =  0.3275911;
		// Save the sign of x
		val sign = if (x < 0) -1 else 1
		val absx =  math.abs(x)

		// A&S formula 7.1.26, rational approximation of error function
		val t = 1.0/(1.0 + p*absx);
		val y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x);
		sign*y
	}

	//from E2LSH code in SelfTuring.cpp
	//w: bin width
	//c: | vector_x - vector_y | , i.e. distance of 2 point
	//it is OK
	def computeFunctionP(w: Double, c: Double):Double = {
		  val x = w / c;
			val M_SQRT2 = math.sqrt(2.0)
			val M_2_SQRTPI = 2.0 / math.sqrt(math.Pi)
			ERF(x / M_SQRT2) - M_2_SQRTPI / M_SQRT2 / x * (1 - math.exp(-(x*x) / 2));
	}

	//from E2LSH & manual
	def computeLfromKP(k: Int, successProb: Double):Int = {
			val p1 = computeFunctionP(binWidthW, radiusR)
		  math.ceil(math.log(1 - successProb) / math.log(1 - math.pow(p1, k))).toInt;
	}

	//from the original paper
	//   & http://cseweb.ucsd.edu/~dasgupta/254-embeddings/lawrence.pdf
	def computeK(p2: Double,n: Int): Int = {
		(math.log(n)/math.log(1.0/p2)).toInt
	}


	def main(args: Array[String]){
		//val baseVector = Source.fromFile("../dataset/sift/siftsmall/test_base").getLines.toList.map{_.split(" ").toIndexedSeq}
		val queryVectors = Source.fromFile("../dataset/sift/siftsmall/query.input").getLines.toList.map{_.split(" ").toIndexedSeq}
		val learnVectors = Source.fromFile("../dataset/sift/siftsmall/learn.input").getLines.toList.map{_.split(" ").map{_.toDouble}.toIndexedSeq}

		val numberOfVectors = learnVectors.length
		val dimension = learnVectors(1).length

		//from the original paper & we let p2 = p(1.0)
		val dimensionK = computeK(computeFunctionP(binWidthW, radiusR), numberOfVectors)
		val numberOfTablesL = computeLfromKP(dimensionK,successProbability)
		val tableSize = numberOfVectors

		//?The value will be minus because radius is less than bin width
		//println("p(R=1.0) = computeFunctionP(4.0,1.0) = " + computeFunctionP(4.0,1.0))


		//learnVectors: 25000 * 128 List

		val hashFunctionH1 = IndexedSeq.fill(dimensionK)((Random.nextDouble()*prime).toLong)
		val hashFunctionH2 = IndexedSeq.fill(dimensionK)((Random.nextDouble()*prime).toLong)
		
		var a = 0
		
		for (vector <- learnVectors){

			def lshFunctionH(): Int = {
				val vectorA = IndexedSeq.fill(dimension)(Random.nextGaussian())
				
				val parameterB = Random.nextDouble() * binWidthW

				val dot = (vector zip vectorA).map{case (a,b) => a * b}.foldLeft(0.0)(_+_)

				math.floor((dot + parameterB)/binWidthW).toInt
			}

			val vectorG = List.fill(dimensionK)(lshFunctionH)

			def computePrimaryHash(h1: IndexedSeq[Long]) = {
				val dot = (vectorG zip h1).map{case (a,b) => a * b}.foldLeft(0L)(_+_)
				val hashValue: Int = ((dot % prime)%tableSize).toInt
				println("hash value: " + hashValue.getClass +" " +hashValue)
			}

			def computeSecondaryHash(h2: IndexedSeq[Long]) = {
				val dot = (vectorG zip h2).map{case (a,b) => a * b}.foldLeft(0L)(_+_)
				val hashValue: Int = (dot % prime).toInt
				println("hash value: " + hashValue.getClass +" "+hashValue)
			}

			if (a < 10){
				computePrimaryHash(hashFunctionH1)
				computeSecondaryHash(hashFunctionH2)
				println
				a += 1
			}
				
		}
		

	}
}
