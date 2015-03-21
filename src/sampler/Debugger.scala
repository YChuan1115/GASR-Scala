package sampler

/**
 * Created by Yang Song on 2015/2/5.
 *
 * This class was created for convenience of debugging
 */
class Debugger(from: Int, to: Int, m: Int, n: Int, r: Int, a: Double = 100.0, b: Double = 100.0)
  extends Sampler(from, to, m, n, r, a, b) {
  def RealInfo() = {
    var sum = 0.0
    var N = 0
    for (item <- real; i = item._1._1; j = item._1._2; rate = item._2) {
      N += 1
      sum += rate * rate
    }
    println(s"Squared sum: $sum Number of Items: $N")
  }
}
