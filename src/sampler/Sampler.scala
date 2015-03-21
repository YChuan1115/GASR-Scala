package sampler

import breeze.linalg._
import breeze.numerics._
import sampler.Distributions._

import scala.collection.mutable.{Map => mMap}
import scala.io.Source

class Sampler(from: Int, end: Int, m: Int, n: Int, r: Int, var a: Double = 100.0, var b: Double = 100.0) {
  //The observed matrix
  val X = Array.ofDim[Double](m + 1, n + 1)
  //The predicted completed matrix and matrix needed to compute averages
  var Z = Array.ofDim[Double](m + 1, n + 1)
  val Sumz = Array.ofDim[Double](m + 1, n + 1)
  //The initial singular vector matrices
  val initlen: Double = 0.9
  // the length of the initial singular vectors
  val U = unitUniform(m, r, initlen)
  val V = unitUniform(n, r, initlen)
  //The squared norms of each singular vector
  val sumU = Array.fill[Double](r + 1)(initlen * initlen)
  val sumV = Array.fill[Double](r + 1)(initlen * initlen)
  //The singular value vectors
  val d = Array.ofDim[Double](r + 1)
  val Sumd = Array.ofDim[Double](r + 1)
  //The gamma values
  val gamma = Array.ofDim[Double](r + 1)
  // The lambda values
  var lambda: Double = 1.0

  var alpha = 0.0
  var beta: Double = 0.0
  //The maps to store the Omega_ij
  val row: mMap[Int, Set[Int]] = mMap.empty
  val column: mMap[Int, Set[Int]] = mMap.empty
  var omega: List[(Int, Int)] = Nil
  var omegaT: List[(Int,Int)] = Nil
  var observations: Int = _
  //The maps storing the true results
  val real: mMap[(Int, Int), Double] = mMap.empty

  def read(files: Array[String], test: Boolean = false): Unit = {
    //Read training or testing data
    files.foreach {
      name => {
        test match {
          case false =>
            Source.fromFile(name).getLines().foreach {
              line => {
                val s = line.split("\t")
                val user = s(0).toInt
                val
                movie = s(1).toInt
                val rate = s(2).toDouble
                X(user)(movie) = rate
                omega = (user -> movie) :: omega
                row.get(user) match {
                  case None => row += user -> Set(movie)
                  case Some(set) => row(user) = set + movie
                }
                column.get(movie) match {
                  case None => column += movie -> Set(user)
                  case Some(set) => column(movie) = set + user
                }
              }
            }
          case true =>
            Source.fromFile(name).getLines().foreach {
              line => {
                val s = line.split("\t")
                val user = s(0).toInt
                val movie = s(1).toInt
                val rate = s(2).toDouble
                val tuple = user -> movie
                real.get(tuple) match {
                  case None =>
                    real += tuple -> rate
                    omegaT = (user -> movie) :: omegaT
                  case _ =>
                    throw new Exception( """
                          Warning! Multiple test key-value pairs!""")
                }
              }
            }
        }
      }
    }
    observations = omega.size
  }


  def init(iter: Int):Unit = {
    /*
      Use MAP algorithm for initialization.
    */
    val w = Array.fill[Double](r + 1)(0.0)
    val flag = Array.fill[Boolean](m + 1, n + 1)(false)

    for ((i, j) <- omega) flag(i - 1)(j - 1) = true

    val matX = DenseMatrix.tabulate(m, n)((i, j) => X(i + 1)(j + 1))
    var matZ = DenseMatrix.tabulate(m, n)((i, j) => Z(i + 1)(j + 1))

    lambda = 1.0

    for (c <- 1 to iter) {
      (1 to r).foreach(t => w(t) = (a + 1) / (b + d(t)))
      for (i <- 0 until m; j <- 0 until n)
        if (flag(i)(j)) matZ(i, j) = matX(i, j)

      val svd.SVD(u, s, v) = svd(matZ)

      for (i <- 0 until r) s(i) = max(0.0, s(i) - 1 / lambda * w(i + 1))

      (1 to r).foreach(i => d(i) = s(i - 1))

      matZ = u(::, 0 until r) * diag(s(0 until r)) * v(0 until r,::)
    }
    val svd.SVD(u, s, v) = svd(matZ)

    for(i <- 0 until m; k <- 0 until r) U(i + 1)(k + 1) = u(i,k) * initlen
    for(i <- 0 until n; k <- 0 until r) V(i + 1)(k + 1) = v(k,i) * initlen
    (0 until r).foreach(i => d(i + 1) = s(i))
    //
    //println(s"Debug: The initial length!: ${sqrt((1 to m).foldLeft(0.0)((s,i)=>s+U(i)(1)*U(i)(1)))}")
    //
  }

  def sampleLambda(): Unit = {
    var sum: Double = 0.0
    for ((i, j) <- omega) {
      sum += (X(i)(j) - Z(i)(j)) * (X(i)(j) - Z(i)(j))
    }
    sum /= 2
    lambda = gamrnd(alpha + observations / 2.0, sum + beta)
  }

  def sampleGammas(): Unit = {
    for (i <- 1 to r)
      gamma(i) = gamrnd(a + 1, b + d(i))
  }

  def sampleds(): Unit = {
    for (a <- 1 to r) {
      val A = omega.foldLeft(0.0) {
        (s, item) => s + U(item._1)(a) * U(item._1)(a) * V(item._2)(a) * V(item._2)(a)
      }

      val B = omega.foldLeft(0.0) {
        (s, item) => {
          var sum = 0.0
          for (k <- 1 to r; if k != a) sum += d(k) * U(item._1)(a) * U(item._1)(k) * V(item._2)(a) * V(item._2)(k)
          s + sum - X(item._1)(item._2) * U(item._1)(a) * V(item._2)(a)
        }
      } + gamma(a) / lambda

      if (-B / A < 0)
        d(a) = tailNorm(mu = -B / A, variance = 1 / A / lambda)
      else {
        d(a) = halfNorm(mu = -B / A, variance = 1 / A / lambda)
      }
    }
  }

  def sampleU(): Unit = {
    /*
      !!!!!!!!!!!!!!!!!!
      Mind the danger that some b never exists in the training data!
      Then U obey the prior distribution, e.g. truncated uniform distribution.
    */
    for (a <- 1 to r; b <- 1 to m) {
      sumU(a) -= U(b)(a) * U(b)(a)
      val bound = sqrt(1 - sumU(a))
      row.get(b) match {
        case Some(set) =>
          val C = set.foldLeft(0.0) {
            (s, j) => s + d(a) * d(a) * V(j)(a) * V(j)(a)
          }
          val D = set.foldLeft(0.0) {
            (s, j) => {
              var sum = 0.0
              for (k <- 1 to r
                   if k != a) sum += d(a) * d(k) * U(b)(k) * V(j)(a) * V(j)(k)
              s + sum - d(a) * X(b)(j) * V(j)(a)
            }
          }
          U(b)(a) = tnorm(mu = -D / C, variance = 1 / C / lambda, a = -bound, b = bound)
        case None =>
          U(b)(a) = tuniform(-bound, bound)
      }
      sumU(a) += U(b)(a) * U(b)(a)
    }

  }

  def sampleV(): Unit = {
    /*
      !!!!!!!!!!!!!!!!!
      Mind the danger that some b never exists in the training data!
      Then V obey the prior distribution, e.g. truncated uniform distribution.
    */

    for (a <- 1 to r; b <- 1 to n) {
      sumV(a) -= V(b)(a) * V(b)(a)
      val bound = sqrt(1 - sumV(a))
      column.get(b) match {
        case Some(set) =>
          val E = set.foldLeft(0.0) {
            (s, i) => s + d(a) * d(a) * U(i)(a) * U(i)(a)
          }
          val F = set.foldLeft(0.0) {
            (s, i) => {
              var sum = 0.0
              for (k <- 1 to r; if k != a) sum += d(a) * d(k) * U(i)(a) * U(i)(k) * V(b)(k)
              s + sum - d(a) * X(i)(b) * U(i)(a)
            }
          }
          V(b)(a) = tnorm(mu = -F / E, variance = 1 / E / lambda, a = -bound, b = bound)
        case None => V(b)(a) = tuniform(-bound, bound)
      }
      sumV(a) += V(b)(a) * V(b)(a)
    }
  }

  def updateZ(): Unit = {
    //Compute the latent matrix Z
    Z = Array.ofDim[Double](m + 1, n + 1)

    /*This is terribly wrong!!!!!!
    for ((i,j) <- omegaT) {
      for (k <- 1 to r) {
        Z(i)(j) += d(k) * U(i)(k) * V(j)(k)
      }
      if (sampleCounter <= end && sampleCounter >= from) {
        Sumz(i)(j) += Z(i)(j)
      }

    }
    if (sampleCounter <= end && sampleCounter >= from) {
      (1 to r).foreach(i => Sumd(i) += d(i))
    }*/
    /* This is also terribly wrong! mind repetitions.
    val needed = omegaT ++ omega
    for ((i,j) <- needed) {
      for (k <- 1 to r) {
        Z(i)(j) += d(k) * U(i)(k) * V(j)(k)
      }
      if (sampleCounter <= end && sampleCounter >= from) {
        Sumz(i)(j) += Z(i)(j)
      }

    }*/

    for (i <- 1 to m; j <- 1 to n) {
      for (k <- 1 to r) {
        Z(i)(j) += d(k) * U(i)(k) * V(j)(k)
      }
      if (sampleCounter <= end && sampleCounter >= from)
        Sumz(i)(j) += Z(i)(j)
    }
    if (sampleCounter <= end && sampleCounter >= from) {
      (1 to r).foreach(i => Sumd(i) += d(i))
    }
  }

  private def _RMSE(Z: Array[Array[Double]]): Double = {
    //must run predict() first.
    //Calculate the predict accuracy.
    var sum = 0.0
    var N = 0

    for (item <- real; i = item._1._1; j = item._1._2; rate = item._2) {
      val dis = Z(i)(j) - rate
      sum += dis * dis
      N += 1
    }
    sqrt(sum / N)
  }

  def RMSE(): Double = {
    _RMSE(Z)
  }

  def averageRMSE(): Double = {
    (1 to r).foreach(i => Sumd(i) /= (end - from + 1))
    _RMSE(Sumz.map(_.map(_ / (end - from + 1).toDouble)))
  }

  private def _err(Z: Array[Array[Double]]): Double = {
    var sum = 0.0
    var tmp = 0.0

    for (item <- real; i = item._1._1; j = item._1._2; rate = item._2) {
      val dis = Z(i)(j) - rate
      sum += dis * dis
      tmp += rate * rate
    }
    sum / tmp
  }

  def err(): Double = {
    _err(Z)
  }

  def averageErr(): Double = {
    (1 to r).foreach(i => Sumd(i) /= (from - end + 1))
    _err(Sumz.map(_.map(_ / (end - from + 1).toDouble)))
  }

  private def _NMAE(Z: Array[Array[Double]]): Double = {
    var N = 0
    var sum = 0.0
    var maxi = -inf
    var mini = inf
    for (item <- real; i = item._1._1; j = item._1._2; rate = item._2) {
      N += 1
      sum += abs(rate - Z(i)(j))
      maxi = max(rate, maxi)
      mini = min(rate, mini)
    }
    sum / N.toDouble / (maxi - mini)
  }

  def NMAE(): Double = _NMAE(Z)

  def averageNMAE(): Double = {
    (1 to r).foreach(i => Sumd(i) /= (end - from + 1))
    _NMAE(Sumz.map(_.map(_ / (end - from + 1).toDouble)))
  }

  var sampleCounter = 0

  def sample() {
    sampleCounter += 1

    //sampleLambda()
    lambda = 1.0

    sampleGammas()
    sampleds()
    sampleU()
    sampleV()
    updateZ()
  }

  //Debugging only!
  /*
    def loadInit() = {
      //
//      load init values from Matlab,
//      only for debugging purposes.

      val readd = Source.fromFile("E:\\Research\\Urgent\\Work of matrix completion\\MCS\\d.txt")
      val readU = Source.fromFile("E:\\Research\\Urgent\\Work of matrix completion\\MCS\\U.txt")
      val readV = Source.fromFile("E:\\Research\\Urgent\\Work of matrix completion\\MCS\\V.txt")
      readd.getLines().foldLeft(1) {
        (i,line) => {
          d(i) = line.trim.toDouble
          i + 1
        }
      }


      readU.getLines().foldLeft(1) {
        (i,line) => {
          line.trim.split(' ').foldLeft(1) {
            (j, num) => {
              U(i)(j) = num.toDouble
              j + 1
            }
          }
          i + 1
        }
      }


      readV.getLines().foldLeft(1) {
        (i,line) => {
          line.trim.split(' ').foldLeft(1) {
            (j, num) => {
              V(i)(j) = num.toDouble
              j + 1
            }
          }
          i + 1
        }
      }

    }
  */
  //loadInit()
  /////////////
}
