# Draft: issue for pyushkevich/greedy (Cedric posts; written 2026-09-27)

Title: `-rj` Jacobian ignores the image direction matrix (correct with `-jac`)

The Jacobian determinant written by `-rj` in the reslice command is computed on the
composed warp while it is still in physical (LPS, mm) units, so the derivatives are
taken along voxel index axes of vectors expressed in a different frame. The result is
only right when the image direction matrix is the identity in ITK terms (voxel axes
LPS-aligned) and the spacing is 1 mm. The standalone `-jac` command converts the warp
with `PhysicalWarpToVoxelWarp` first and is correct.

Reproduction (greedy 1.4.0.3 via the picsl_greedy wheel, arm64 Linux): a dark sphere of
radius 8 voxels registered to the same sphere at radius 10 (volume ratio 1.95), SSD,
`-n 100x60x30 -e 0.7 -s 2vox 0.5vox`, identical voxel data under three affines:

| voxel-to-world | det from `-rj` inside the sphere | det from `-jac` |
|---|---|---|
| identity (RAS) | 0.675 | 2.033 |
| diag(-1, 1, 1) (LAS) | 1.183 | not run |
| diag(-1, -1, 1) (LPS-aligned) | 2.034 | not run |

With a radial stretch 1+k, `-rj` gives (1-k)^2 (1+k) for RAS and (1-k)(1+k)^2 for LAS in
place of (1+k)^3, i.e. the sign flips per flipped axis; NIfTI files with RAS or LAS
affines (most nibabel/FSL output) get a wrong, sometimes inverted, Jacobian.

Where: `GreedyAPI.cxx`, `RunReslice`, block `if(r_param.out_jacobian_image.size())`
calls `LDDMMType::field_jacobian_det(warp, iTemp)` on the physical warp returned by
`ReadTransformChain`. `RunJacobian` does `OFHelperType::PhysicalWarpToVoxelWarp(warp, warp, warp)`
before differentiating.

Suggested fix (one extra conversion on a copy, the physical warp is still needed for
the reslicing that follows):

```cpp
  if(r_param.out_jacobian_image.size())
    {
      VectorImagePointer warp_vox = LDDMMType::new_vimg(warp);
      OFHelperType::PhysicalWarpToVoxelWarp(warp, warp, warp_vox);
      ImagePointer iTemp = ImageType::New();
      LDDMMType::alloc_img(iTemp, warp);
      LDDMMType::field_jacobian_det(warp_vox, iTemp);
      WriteImageViaCache(iTemp.GetPointer(), r_param.out_jacobian_image.c_str(), itk::IOComponentEnum::FLOAT);
    }
```

Workaround for users: write the composed warp with `-rc` and run `greedy -d 3 -jac warp out`.

Found while building lesiontrack (MS lesion Jacobian tracking); its own determinant is
now computed in the voxel frame and matches `-jac` to three decimals on the phantom.
