{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE BangPatterns #-}

module Lib
    ( trilinear
    , resample
    , affine
    , array_correlation
    , fisher_transform
    , naive_12p_registration
    ) where

import Data.Array.Repa
import qualified Data.Array.Repa as Repa
import Data.Array.Repa.Algorithms.Matrix
import Statistics.Sample
import Data.Array.Repa.Repr.Unboxed
import Data.Array.Repa.Repr.Vector

type DIM3D = (Z :. Double :. Double :. Double)

trilinear :: (Floating a, Fractional a) => Array V DIM3 a -> DIM3D -> a
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
        ((fromIntegral x0) - (realToFrac x)) / (fromIntegral (x0 - x1)) else
        0
      yd = if y0 /= y1 then
        ((fromIntegral y0) - (realToFrac y)) / (fromIntegral (y0 - y1)) else
        0
      zd = if z0 /= z1 then
        ((fromIntegral z0) - realToFrac z) / (fromIntegral (z0 - z1)) else
        0
      in let
        c00 = v000 * (1 - xd) + v100 * xd
        c01 = v001 * (1 - xd) + v101 * xd
        c10 = v010 * (1 - xd) + v110 * xd
        c11 = v011 * (1 - xd) + v111 * xd
        in let
          c0 = c00 * (1 - yd) + c10 * yd
          c1 = c01 * (1 - yd) + c11 * yd
          in  c0 * (1 - zd) + c1 * zd


resample ::
  (Floating a, Floating b, Fractional a, Fractional b) =>
  Array V DIM3 a ->
  (DIM3 -> DIM3D) -> --fun from coordinate system to coordinate system
  (Array V DIM3 a -> DIM3D -> b) -> --interpolator
  Array D DIM3 b
  
resample !arr !xfm !interp =
  Repa.traverse arr id xfm_and_interp
  where
    xfm_and_interp = \_ d -> ((interp arr (xfm d)))


affine ::
  (Floating a, Fractional a, Real a) =>
  DIM3 -> Array V DIM2 a -> DIM3D

affine !(Z :. x :. y  :. z) !xfm =
  (Z :. x' :. y' :. z')
  where
    [x',y',z',_] = (toList (mmultS (computeS $ Repa.map realToFrac xfm) pnt_array))
      where
        pnt_array = (fromListUnboxed
                     (Z :. (3 :: Int) :. (1 :: Int))
                     (fmap fromIntegral [x,y,z])) -- :: Array U DIM2 Double


array_correlation ::
  (Real a, Fractional a) => Array V DIM3 a -> Array V DIM3 a -> a
array_correlation !arr1 !arr2 =
  realToFrac
   (correlation
    (toUnboxed
     (computeS
      (Repa.zipWith coerce_and_tuple
         arr1 arr2))))
  where
     coerce_and_tuple x y = (realToFrac x, realToFrac y) 


fisher_transform = atanh

-- Test code
naive_12p_registration ::
  (Floating a, Fractional a, Real a) =>
  Array V DIM3 a -> [a] -> a

naive_12p_registration !arr !params =
  array_correlation arr (computeS $ resample arr (flip affine xfm) trilinear)
  where
    xfm = fromListVector (Z :. 4 :. 4) params
