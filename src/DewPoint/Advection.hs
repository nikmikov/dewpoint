module DewPoint.Advection where



-- | parameter advection due to winds
--
--     -----> U
--
--       dx
--  V  -------
--     |_____|_______  _
--  ^  |     |       |  |
--  |  | cXY |  cX   |  |  dy
--  |  |_____|_______| _|
--     |     |       |
--     | cY  |  c0   |
--     |     |       |
--     ---------------
--
--
-- return pair of coefficients indicating  cell quantity that will be replaced
--  with content from source cells by (X, Y) axis
--  values are in range 0..1
advectionCoefficients :: (Floating a, Ord a) =>
                      a -- ^ U-wind along X (m/s)
                   -> a -- ^ V-wind along Y (m/s)
                   -> a -- ^ length of cell side (assuming cell is always a square)
                   -> a -- ^ time interval (sec)
                   -> (a, a)
advectionCoefficients _ _ 0 _ = (0, 0)
advectionCoefficients u v l dt =
    let dx = min (abs(u) * dt) l
        dy = min (abs(v) * dt) l
        totalArea = l * l
        cXY = dx * dy / totalArea
        cX = (l - dy) * dx / totalArea
        cY = (l - dx) * dy / totalArea
        cX' = if (dx + dy > 0) then dx / (dx + dy) else 0
        cY' = if (dx + dy > 0) then 1 - cX' else 0
        rX = cX + cX' * cXY
        rY = cY + cY' * cXY
        rY' = if rX + rY > 1.0 then 1.0 - rX else rY
    in (rX, rY')
