module Main where

import Lib
import Data.Array.Repa
import qualified Data.Array.Repa as Repa

main :: IO ()
main =
  do let x = fromListUnboxed (Z :. (3::Int) :. (3::Int) :. (3::Int)) [1..27] ::
           Array U DIM3 Double
     let xfm = fromListUnboxed
              (Z :. (4::Int) :. (4::Int))
              [0.9, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]    
     let res = computeS $ resample x (flip affine xfm) trilinear :: Array U DIM3 Double
     print res
     
