package potrace;


import org.opencv.core.Mat;
import org.opencv.imgcodecs.Imgcodecs;
import org.opencv.imgproc.Imgproc;

import java.awt.*;
import java.util.ArrayList;


public class Potrace {


    /**
     * This function converts source image into a boolean array
     *
     * @param srcImage
     * @return the boolean array where true stands for background and false
     * stands for foreground
     * @apiNote the output array has larger size than source image cause
     * background border is added around the image by 1 pixel
     */
    public static boolean[][] bitmapToBinary(Mat srcImage) {

        assert (srcImage != null && !srcImage.empty());

        Mat binaryImage = new Mat();
        if (srcImage.channels() == 3) {
            Mat grayImage = new Mat();
            Imgproc.cvtColor(srcImage, grayImage, Imgproc.COLOR_BGR2GRAY);
            Imgproc.threshold(grayImage, binaryImage,
                    207, 255, Imgproc.THRESH_BINARY);
        } else {
            Imgproc.threshold(srcImage, binaryImage,
                    207, 255, Imgproc.THRESH_BINARY);
        }

        boolean[][] result = new boolean[srcImage.rows() + 2][srcImage.cols() + 2];
        for (int row = 0; row < binaryImage.rows(); row++) {
            for (int col = 0; col < binaryImage.cols(); col++) {
                if ((int) binaryImage.get(row, col)[0] == 255) {
                    result[row + 1][col + 1] = true;
                } else {
                    result[row + 1][col + 1] = false;
                }
            }
        }

        // add white borders
        for (int col = 0; col < result[0].length; col++) {
            result[0][col] = true;
            result[result.length - 1][col] = true;
        }
        for (int row = 1; row < result.length - 1; row++) {
            result[row][0] = true;
            result[row][result[0].length - 1] = true;
        }

        return result;
    }


    /**
     * This function attempts to trace the given bitmap using the
     * given tracing parameters.
     *
     * @param bm                the given bitmap in binary array form
     * @param ListOfCurveArrays the curve result of Potrace
     */
    public static void potrace_trace(boolean[][] bm, ArrayList ListOfCurveArrays) {

        //Step 1
        ArrayList plistp = new ArrayList();
        bm_to_pathlist(bm, plistp);

        //Step 2
        process_path(plistp);

        //Step 3
        PathList_to_ListOfCurveArrays(plistp, ListOfCurveArrays);
    }


    /**
     * Decompose the given bitmap into paths. Returns a list of Path
     * objects with the fields len, pt, area filled
     *
     * @param bm     A binary bitmap which holds the imageinformations.
     * @param plistp List of Path objects
     */
    private static void bm_to_pathlist(boolean[][] bm, ArrayList plistp) {

        Integer x = 0, y = 0;
        Point temp;

        while (true) {
            temp = findNext(bm, x, y);
            if (temp == null) {
                break;
            } else {
                getContour(bm, temp.x, temp.y, plistp);
            }
        }

    }

    private static void getContour(boolean[][] bm, Integer x, Integer y, ArrayList plistp) {

        Path Contur = findpath( bm, new IntPoint( x, y ) );

        Xor_Path( bm, Contur );
        ArrayList PolyPath = new ArrayList();
        // only area > turdsize is taken

        if( Contur.area > turdsize )
        {
            plistp.Add( PolyPath );
            PolyPath.Add( Contur ); // Path with index 0 is a conture
        }


        while( FindNext( bm, ref x, ref y, Contur ) )
        {
            Path Hole = findpath( bm, new IntPoint( x, y ) );
            //Path Hole = findpath(bm, x, y);
            Xor_Path( bm, Hole );
            if( Hole.area > turdsize )

                PolyPath.Add( Hole ); // Path with index > 0 is a hole,
            while( FindNext( bm, ref x, ref y, Hole ) ) // 13.07.12 von if auf while
                getContur( bm, x, y, plistp );

        }
    }


    /**
     * Searches a x and a y such that source[x,y] = true and source[x+1,y] false.
     * If this not exists, false will be returned else the result is true.
     *
     * @param bm a binary matrix
     * @param x  x index in the source Matrix
     * @param y  y index in the source Matrix
     * @return
     */
    private static Point findNext(boolean[][] bm, Integer x, Integer y) {

        for (int col = 1; col < bm.length - 1; col++) {
            for (int row = 0; row < bm.length - 1; row++) {
                // black found
                if (!bm[row + 1][col]) {
                    x = row;
                    y = col;
                    return new Point(row, col);
                }
            }
        }

        return null;

    }


    private static void process_path(ArrayList plistp) {

    }

    private static void PathList_to_ListOfCurveArrays(ArrayList plistp, ArrayList ListOfCurveArrays) {

    }

    //Below is for test
    public static void main(String[] args) {

        String dllPath = "C:\\OpenCV\\opencv\\build\\java\\x64\\opencv_java320.dll";
        System.load(dllPath);

        testBitmapToBinary();
    }

    private static void testBitmapToBinary() {

        String filePath = "E:\\Java_Projects\\Potrace\\resources\\sourceEntireImages\\4a.png";
        Mat srcImage = Imgcodecs.imread(filePath);
        boolean[][] result = bitmapToBinary(srcImage);

        for (int i = 0; i < result.length; i++) {
            for (int j = 0; j < result[i].length; j++) {
                if (result[i][j]) {
                    System.out.print(0);
                } else {
                    System.out.print(1);
                }
            }
            System.out.println();
        }
    }
}
