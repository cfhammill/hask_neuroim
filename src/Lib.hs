{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE BangPatterns #-}

module Lib
    ( trilinear
    , resample
    , affine
    ) where

import Data.Array.Repa
import qualified Data.Array.Repa as Repa
import Data.Array.Repa.Algorithms.Matrix
import Statistics.Sample

type DIM3D = (Z :. Double :. Double :. Double)

trilinear :: Array U DIM3 Double -> DIM3D -> Double
trilinear !arr !(Z :. x :. y :. z) =
  let
    x0 = floor x
    x1 = ceiling x :: Int
    y0 = floor y
    y1 = ceiling y :: Int
    z0 = floor z
    z1 = ceiling z :: Int
    in let
      v000 = arr ! (Z :. x0 :. y0 :. z0)
      v100 = arr ! (Z :. x1 :. y0 :. z0)
      v010 = arr ! (Z :. x0 :. y1 :. z0)
      v110 = arr ! (Z :. x1 :. y1 :. z0)
      v001 = arr ! (Z :. x0 :. y0 :. z1)
      v101 = arr ! (Z :. x1 :. y0 :. y1)
      v011 = arr ! (Z :. x0 :. y1 :. y1)
      v111 = arr ! (Z :. x1 :. y1 :. y1)
      xd = if x0 /= x1 then
        ((fromIntegral x0) - x) / (fromIntegral (x0 - x1)) else
        0
      yd = if y0 /= y1 then
        ((fromIntegral y0) - y) / (fromIntegral (y0 - y1)) else
        0
      zd = if z0 /= z1 then
        ((fromIntegral z0) - z) / (fromIntegral (z0 - z1)) else
        0
      in let
        c00 = v000 * (1 - xd) + v100 * xd
        c01 = v001 * (1 - xd) + v101 * xd
        c10 = v010 * (1 - xd) + v110 * xd
        c11 = v011 * (1 - xd) + v111 * xd
        in let
          c0 = c00 * (1 - yd) + c10 * yd
          c1 = c01 * (1 - yd) + c11 * yd
          in   c0 * (1 - zd) + c1 * zd


resample ::
  Array U DIM3 Double ->
  (DIM3 -> DIM3D) -> --fun from coordinate system to coordinate system
  (Array U DIM3 Double -> DIM3D -> Double) -> --interpolator
  Array D DIM3 Double
  
resample !arr !xfm !interp =
  Repa.traverse arr id xfm_and_interp
  where
    xfm_and_interp = \_ d -> ((interp arr (xfm d)))


affine ::
  DIM3 -> Array U DIM2 Double -> DIM3D

affine !(Z :. x :. y  :. z) !xfm =
  (Z :. x' :. y' :. z')
  where
    [x',y',z',_] = (toList (mmultS xfm pnt_array))
      where
        pnt_array = (fromListUnboxed
                     (Z :. (3 :: Int) :. (1 :: Int))
                     (fmap fromIntegral [x,y,z])) -- :: Array U DIM2 Double


array_correlation :: Array U DIM3 Double -> Array U DIM3 Double -> Double
array_correlation !arr1 !arr2 =
  correlation (toUnboxed (computeS (Repa.zipWith (,) arr1 arr2)))

fisher_transform :: Double -> Double
fisher_transform x = 0.5 * log ( (1 + x) / (1 - x) )
