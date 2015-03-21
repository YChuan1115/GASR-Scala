package sampler

import breeze.numerics.{digamma, abs, log, trigamma}
import Distributions._

/**
 * Created by Yang Song on 2015/2/4.
 *
 * Department of Physics, Tsinghua University
 * Beijing, China
 *
 */
class EMSampler(from: Int, end: Int, m: Int, n: Int,
                r: Int, _a: Double = 100.0,
                _b: Double = 100.0, EMNumber: Int = 5) extends Sampler(from, end, m, n, r, _a, _b) {
  //Scala do not allow supper.a or supper.b temporarily, extremely crazy! {
  val S = Array.ofDim[Double](EMNumber + 1)
  val SL = Array.ofDim[Double](EMNumber + 1)
  val SLLambda = Array.fill(EMNumber + 1)(0.0)
  val SLambda = Array.fill(EMNumber + 1)(0.0)
  var EMCounter = 0
  //Note the values to remove

  override def sampleGammas(): Unit = {
    for (i <- 1 to r) {
      gamma(i) = gamrnd(a + 1, b + d(i))
    }
    //start up stochastic EM!
    S(EMCounter) = gamma.drop(1).sum
    SL(EMCounter) = gamma.drop(1).map(log(_)).sum
  }

  override def sampleLambda(): Unit = {
    var sum: Double = 0.0
    for ((i, j) <- omega) {
      sum += (X(i)(j) - Z(i)(j)) * (X(i)(j) - Z(i)(j))
    }
    sum /= 2
    lambda = gamrnd(alpha + observations / 2.0, sum + beta)

    println(s"Effective alpha: ${alpha+observations/2.0}, effective beta:${sum+beta}")
    SLambda(EMCounter) = lambda
    SLLambda(EMCounter) = log(lambda)
  }

  private def EM(): Unit = {

    val s = S.drop(1).sum / EMNumber
    val sl = SL.drop(1).sum / EMNumber
    val slambda = SLambda.drop(1).sum / EMNumber
    val sllambda = SLLambda.drop(1).sum / EMNumber
    var ta = a + 1.0
    //The initial values has been changed compared with the original Matlab codes.
    //Deprecated: a = 0.1

    while (abs(ta - a) > 1e-3) {
      ta = a
      a -= (digamma(a) - log(r * a / s) - sl / r) / (trigamma(a) - 1 / a)
      if (a < 0) a = 0.1
    }
    b = r * a / s

    alpha = 0.1
    var talpha = alpha + 1.0
    while(abs(talpha - alpha) > 1e-3){
      talpha = alpha
      alpha -= (log(alpha) - digamma(alpha) + sllambda - log(slambda))/(1.toDouble/alpha - trigamma(alpha))
      if(alpha < 0) alpha = 0.1
      if(alpha >= 400000) alpha = 400000
    }
    beta = alpha / slambda
  }

  override def sample(): Unit = {

    EMCounter += 1
    sampleCounter += 1

    sampleLambda()
    sampleGammas()
    sampleds()
    sampleU()
    sampleV()
    updateZ()

    if (EMCounter >= EMNumber) {
      EM()
      EMCounter = 0
    }
  }
}
