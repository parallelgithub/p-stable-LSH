import scala.io.Source
import scala.util.Random

object LSH {
	def main(args: Array[String]){
		//val baseVector = Source.fromFile("../dataset/sift/siftsmall/test_base").getLines.toList.map{_.split(" ").toIndexedSeq}
		val queryVectors = Source.fromFile("../dataset/sift/siftsmall/query.input").getLines.toList.map{_.split(" ").toIndexedSeq}
		val learnVectors = Source.fromFile("../dataset/sift/siftsmall/learn.input").getLines.toList.map{_.split(" ").map{_.toDouble}.toIndexedSeq}

		val dimention = learnVectors(1).length

		//learnVectors: 25000 * 128 List
		/*
		for (query <- queryVectors){
			if(learnVectors.contains(query))
					println("Query #"+queryVectors.indexOf(query)+" is in learnVector #"+learnVectors.indexOf(query))
		}
		*/
		
		val parameterW = 4.0
		var a = 0
		for (vector <- learnVectors){

			if(a<20){
			val vectorA = IndexedSeq.fill(dimention)(Random.nextGaussian())
			
			val parameterB = Random.nextDouble() * parameterW

			val dot = (vector zip vectorA).map{case (a,b) => a * b}.foldLeft(0.0)(_+_)
			val hashValue = ((dot + parameterB)/parameterW).toInt
			println(hashValue)
			

			a += 1
			}
		}

	}
}
