module Main where

import Lib

main :: IO ()
main = do
  let x = fromListUnboxed (Z :. (3::Int) :. (3::Int) :. (3::Int)) [1..27] ::
        Array U DIM3 Double;
  let i4 =
        fromListUnboxed
          (Z :. (4::Int) :. (4::Int))
          [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
  print $ (computeP $ resample x (flip affine i4) trilinear  :: IO (Array U DIM3 Double))
  
