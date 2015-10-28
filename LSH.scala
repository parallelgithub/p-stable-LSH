import scala.io.Source

object LSH {
	def main(args: Array[String]){
		val fileLines = Source.fromFile("../dataset/sift/siftsmall/siftsmall_base.input").getLines.toList.map{_.split(" ").toList}

		//fileLines: 10000 * 128 List

	}
}
