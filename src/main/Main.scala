/**
 * Created by Yang Song on 2015/2/3.
 *
 * Department of Physics
 * Tsinghua University
 */
package main

import java.io.PrintWriter

import breeze.linalg.max
import sampler.{Sampler, EMSampler}

object Main {
  def main(args: Array[String]) {
    //val from = args(0).toInt
    //val end = args(1).toInt
    val N = 10
   // for(N <- from to end) {
    val diary = new PrintWriter(s"diary$N.txt")
    //val diary = new PrintWriter("Movielens100k.txt")
    // MovieLens 100k
/*
              val m = 943
              val n = 1682
              val r = 100
              val iter = 100
              val averageNumber = 40
    /*
              val trainFiles = Array(
                "E:\\Research\\Urgent\\Work of matrix completion\\Movielens\\ml-100k\\ml-100k\\ua.base",
                "E:\\Research\\Urgent\\Work of matrix completion\\Movielens\\ml-100k\\ml-100k\\ub.base")
              val testFiles = Array(
                "E:\\Research\\Urgent\\Work of matrix completion\\Movielens\\ml-100k\\ml-100k\\ua.test",
                "E:\\Research\\Urgent\\Work of matrix completion\\Movielens\\ml-100k\\ml-100k\\ub.test")
                */
    val trainFiles = Array("ua.base")
    val testFiles = Array("ua.test")
*/
      // Complete Case

      val m = 100
      val n = 100
      val r = 100
      val iter = 500
      val averageNumber = 50
      val trainFiles = Array(s"X$N.txt")
      val testFiles = Array(s"Z$N.txt")


          // MovieLens 1M
   /*
          val m = 6040
          val n = 3952
          val r = 30
          val iter = 100
          val averageNumber = 40
          */
      //    val trainFiles = Array("E:\\Research\\Urgent\\Work of matrix completion\\Movielens\\ml-1m\\ml-1m\\train.dat")
       //   val testFiles = Array("E:\\Research\\Urgent\\Work of matrix completion\\Movielens\\ml-1m\\ml-1m\\test.dat")
      //    val trainFiles = Array("train.dat")
       //   val testFiles = Array("test.dat")

      /*
              // Jester 1
              val m = 24983
              val n = 100
              val r = 40
              val iter = 100
              val averageNumber = 50
              val trainFiles = Array("train1.dat")
              val testFiles = Array("test1.dat")
      */

      // Jester 2
      /*
          val m = 23500
          val n = 100
          val r = 40
          val iter = 100
          val averageNumber = 50
          val trainFiles = Array("train2.dat")
          val testFiles = Array("test2.dat")
      */
      /*
          // Jester 3
          val m = 24938
          val n = 100
          val r = 40
          val iter = 100
          val averageNumber = 50
          val trainFiles = Array("train3.dat")
          val testFiles = Array("test3.dat")
      */
      val learner = new EMSampler(max(iter - averageNumber + 1, 0), iter, m = m, n = n, r = r,EMNumber = 3)
      learner.read(trainFiles.take(1))
      learner.read(testFiles.take(1), test = true)
      //learner.init(1)

      val t1 = System.nanoTime()

      for (i <- 1 to iter) {
        val begin = System.nanoTime()
        learner.sample()
        learner.d.drop(1).sorted.foreach(diary.println)
        learner.d.drop(1).sorted.foreach(println)
        //diary.println(s"a = ${learner.a} b = ${learner.b}")
        println(s"a = ${learner.a} b = ${learner.b}")
        println(s"alpha = ${learner.alpha} beta = ${learner.beta}")
        println(s"lambda = ${learner.lambda}")
        val error = learner.NMAE()
        diary.println(s"Iteration: $i, Accuracy: $error")
        println(s"Iteration: $i, Accuracy: $error")
        println(s"Elasped Time: ${(System.nanoTime() - begin) / 1e9}s")
        diary.println(s"Elasped Time: ${(System.nanoTime() - begin) / 1e9}s")
      }
      val error = learner.averageNMAE()
      println(s"Averaged Accuracy for $averageNumber iterations: $error")
      //println("Averaged d's : ")
      //learner.Sumd.sorted.foreach(println)
      //diary.println("Averaged d's : ")
      learner.Sumd.drop(1).sorted.foreach(diary.println)
      diary.println(s"Averaged Accuracy for $averageNumber iterations: $error")
      diary.println(s"Averaged Elapsed Time: ${(System.nanoTime() - t1) / iter / 1e9}s")
      diary.close()
    //}
  }
}