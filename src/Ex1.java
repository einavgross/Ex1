/**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
 * https://docs.google.com/document/d/1GcNQht9rsVVSt153Y8pFPqXJVju56CY4/edit?usp=sharing&ouid=113711744349547563645&rtpof=true&sd=true
 *
 * This class represents a set of static methods on a polynomial functions - represented as an array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynomial function: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code below")
 *
 * @author boaz.benmoshe

 */
public class Ex1 {
	/** Epsilon value for numerical computation, it serves as a "close enough" threshold. */
	public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
	/** The zero polynomial function is represented as an array with a single (0) entry. */
	public static final double[] ZERO = {0};
	/**
	 * Computes the f(x) value of the polynomial function at x.
	 * @param poly - polynomial function
	 * @param x
	 * @return f(x) - the polynomial function value at x.
	 */
	public static double f(double[] poly, double x) {
		double ans = 0;
		for(int i=0;i<poly.length;i++) {
			double c = Math.pow(x, i);
			ans += c*poly[i];
		}
		return ans;
	}
	/** Given a polynomial function (p), a range [x1,x2] and an epsilon eps.
	 * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps, 
	 * assuming p(x1)*p(x2) <= 0.
	 * This function should be implemented recursively.
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
	 */
	public static double root_rec(double[] p, double x1, double x2, double eps) {
		double f1 = f(p,x1);
		double x12 = (x1+x2)/2;
		double f12 = f(p,x12);
		if (Math.abs(f12)<eps) {
            return x12;}
		if(f12*f1<=0) {
            return root_rec(p, x1, x12, eps);}
		else {
            return root_rec(p, x12, x2, eps);}
	}
	/**
	 * This function computes a polynomial representation from a set of 2D points on the polynom.
	 * The solution is based on: //	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
	 * Note: this function only works for a set of points containing up to 3 points, else returns null.
	 * @param xx an array with values of x1,x2,x3
	 * @param yy an array with values of y1,y2,y3
	 * @return an array of doubles representing the coefficients of the polynom.
     * this function uses sub function to find the Polynom equation from 2 or 3 points, based on the formula above
     * double [] ans = null;
     * int lx = xx.length;
     * int ly = yy.length;
     * if(xx!=null && yy!=null && lx==ly && lx>1 && lx<4){ //The arrays must be not null and have the same length (2 or 3)
     *    if(lx==2){ans = polyForTwo(xx,yy);} //use the sub-function for 2 points
     *    if(lx==3){ans=polyForThree(xx,yy);} //use the sub-function for 3 points
     * }
     * return ans //return final array
	 */
	public static double[] PolynomFromPoints(double[] xx, double[] yy) {
        double [] ans = null;
		int lx = xx.length;
		int ly = yy.length;
		if(xx!=null && yy!=null && lx==ly && lx>1 && lx<4) {
            if(lx==2){
                ans = polyForTwo(xx,yy);
            }
            if(lx==3){
                ans=polyForThree(xx,yy);
            }
		}
		return ans;
	}
	/** Two polynomials functions are equal if and only if they have the same values f(x) for n+1 values of x,
	 * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
	 * @param p1 first polynomial function
	 * @param p2 second polynomial function
	 * @return true if p1 represents the same polynomial function as p2.
     * this function creates (with sub-function) an array of n+1 X values.
     * the function return false if one of the values doesn't return the same f(x) value(up to EPS) for p1 and p2.else return true
     * boolean ans=true
     * n=max(p1.length-1,p1.length-1)               //set n to the max of the max index of each array - representing the max degree
     * double [] x = createFinalArr(n+1);           //a sub-function that creates an array with size n+1
     * for (int i = 0; i < x.length && ans; i++) {   //go through the whole array or until false
     *     if (Math.abs(f(p1,x[i])-f(p2,x[i]))>EPS) {ans=false;}    //check if the f(x) for each x is not equals
     * }                                                            //ends for block
     * return ans
	 */
	public static boolean equals(double[] p1, double[] p2) {
		boolean ans = true;
        int n= Math.max(p1.length-1,p2.length-1);
        double [] x = createFinalArr(n+1);
        for (int i = 0; i < x.length && ans; i++) {
            if (Math.abs(f(p1,x[i])-f(p2,x[i]))>EPS){
                ans=false;
            }
        }
		return ans;
	}

	/** 
	 * Computes a String representing the polynomial function.
	 * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
	 * @param poly the polynomial function represented as an array of doubles
	 * @return String representing the polynomial function:
     * this function concats each value in the array (using StringBuilder()) with the right sign (+/-) and the power (x^Index)
     * the function skips values with 0 in the array and ignores ^1 (x).
     *StringBuilder ans=newStringBuilder();
     * if(poly.length==0){ans=newStringBuilder("0");}
     * 	else{
     *   	String val="";
     *   	for(inti=poly.length-1;i>=0;i--){
     * 	        if(poly[i]!=0){                                             //append only coefficients !=0
     * 		        if(poly[i]<0|val.isEmpty()) {val=poly[i]}               //append with no sign if the coefficient is negative or first
     * 		        if(poly[i]>0&&!val.equals(poly[i])) {val= "+"+poly[i]}  //append + if the coefficient is positive
     * 		        if(i>1){ans.append(val).append("x^").append(i)}         //append the power to x^ if not 0 or 1
     * 		        else{                                                   //if i=1 or i=o
     * 		            if(i==1){ans.append(val).append("x");}              //if the power is 1, append the coefficient to x
     * 		            if(i==0) {ans.append(val);}                         //if the power is 0, append the coefficient alone
     *              }                                                       //ends else block
     *          }                                                           //ends if block
     *      }                                                               //ends for block
     * }                                                                    //ends else block
     * return ans.toString();                                               //convert ans to regular string and return
	 */
	public static String poly(double[] poly) {
		StringBuilder ans = new StringBuilder();
		if(poly.length==0) {
            ans = new StringBuilder("0");}
		else {
            String val = "";
            for(int i=poly.length-1;i>=0;i--) {
                if (poly[i]!=0) {
                    if (poly[i] < 0 | val.isEmpty()) {
                        val = ""+poly[i];
                    }
                    if (poly[i] > 0 && !val.equals(""+poly[i])) {
                        val = " +" + poly[i];
                    }
                    if (i>1) {
                        ans.append(val).append("x^").append(i);
                    }
                    else {
                        if (i==1){
                            ans.append(val).append("x");
                        }
                        if (i==0){
                            ans.append(val);
                        }
                    }
                }
            }
		}
		return ans.toString();
	}
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an epsilon eps. This function computes an x value (x1<=x<=x2)
	 * for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
	 */
	public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
		double ans = x1;
        /** add you code below

         /////////////////// */
		return ans;
	}
	/**
	 * Given a polynomial function (p), a range [x1,x2] and an integer with the number (n) of sample points.
	 * This function computes an approximation of the length of the function between f(x1) and f(x2) 
	 * using n inner sample points and computing the segment-path between them.
	 * assuming x1 < x2. 
	 * This function should be implemented iteratively (none recursive).
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfSegments - (A positive integer value (1,2,...).
	 * @return the length approximation of the function between f(x1) and f(x2).
	 */
	public static double length(double[] p, double x1, double x2, int numberOfSegments) {
		double ans = x1;
        /** add you code below

         /////////////////// */
		return ans;
	}
	
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an integer representing the number of Trapezoids between the functions (number of samples in on each polynom).
	 * This function computes an approximation of the area between the polynomial functions within the x-range.
	 * The area is computed using Riemann's like integral (https://en.wikipedia.org/wiki/Riemann_integral)
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfTrapezoid - a natural number representing the number of Trapezoids between x1 and x2.
	 * @return the approximated area between the two polynomial functions within the [x1,x2] range.
	 */
	public static double area(double[] p1,double[]p2, double x1, double x2, int numberOfTrapezoid) {
		double ans = 0;
        /** add you code below

         /////////////////// */
		return ans;
	}
	/**
	 * This function computes the array representation of a polynomial function from a String
	 * representation. Note:given a polynomial function represented as a double array,
	 * getPolynomFromString(poly(p)) should return an array equals to p.
	 * 
	 * @param p - a String representing polynomial function.
	 * @return
	 */
	public static double[] getPolynomFromString(String p) {
		double [] ans = ZERO;//  -1.0x^2 +3.0x +2.0
        /** add you code below

         /////////////////// */
		return ans;
	}
	/**
	 * This function computes the polynomial function which is the sum of two polynomial functions (p1,p2)
	 * @param p1 An array of doubles representing a polynomial function.
	 * @param p2 An array of doubles representing a polynomial function.
	 * @return an array representing the polynomial function of the sum of the two polynomial functions.
     * this function compute the sum of the coefficients for each value and insert it in the new array, by the value's power.
     * double [] ans = ZERO;
     * double [] max_p = maxArray(p1, p2);  // a private function to get the longer array
     * double [] min_p= minArray(p1,p2);    //a private function to get the shorter array
     * ans = new double [max_p.length] //change the ans array size to the length of the longer array
     * for(i=0;i<ans.length;i++){      //go through all the values of the final array
     * if(i<min_p.length) {         //if i hasn't reached the end of the min array => more objects to sum up
     *    ans[i]= max_p[i]+min_p[i];} //sum up the values with the same power (same index)
     * else {                       //no more values to sum in from min array
     *  ans[i] = max_p[i];}       //copy the rest of the values from max array to the final array
     * }                          //  ends for block
     * return ans
	 */
	public static double[] add(double[] p1, double[] p2) {
		double [] ans = ZERO;
        double [] max_p =Ex1.maxArray(p1, p2);
        double [] min_p= Ex1.minArray(p1,p2);
        ans = new double[max_p.length];
        for(int i=0;i<ans.length;i++) {
            if (i<min_p.length) {
                ans[i] = max_p[i]+min_p[i];
            }
            else {
                ans[i] = max_p[i];
            }
        }
		return ans;
	}

    /**
	 * This function computes the polynomial function which is the multiplication of two polynoms (p1,p2)
	 * @param p1 An array of doubles representing a polynomial function.
	 * @param p2 An array of doubles representing a polynomial function.
	 * @return an array representing the polynomial function of the multiplication of the two polynomial functions.
     * this function multiply each value in the first array with each value in the second array and insert it in the new array
     * double [] ans = ZERO;
     * if(!isAllZero(p1) && !isAllZero(p2)) {     //if one of the poly is only zeros, so the answer is also zero
     * ans = new double [(p1.length+p2.length)-1] //set the final array's length to the sum of p1.length and p2.length, minus 1
     * for(int i=0;i<p1.length;i++){               //go through all the values in p1
     *  for(int j=0;j<p2.length;j++){              //go through all the values in p2
     *      ans[i+j]=ans[i+j]+(p1[i]*p2[j])}       //sums up the multiplications of each of the values in the right index
     * }                                            //ends 1st for block
     * }                                            //ends 1st for block
     * return ans
	 */
	public static double[] mul(double[] p1, double[] p2) {
		double [] ans = ZERO;
        if(!isAllZero(p1) && !isAllZero(p2)) {
            ans = new double[(p1.length + p2.length) - 1];
            for (int i = 0; i < p1.length; i++) {
                for (int j = 0; j < p2.length; j++) {
                    ans[i + j] = ans[i + j] + (p1[i] * p2[j]);
                }
            }
        }
		return ans;
	}
	/**
	 * This function computes the derivative of the po polynomial function.
	 * @param po Array of doubles representing the polynomial function.
	 * @return Array of doubles representing the derivative of the polynomial function.
     * The function multiplies the coefficient of X by its power and inserts it into the correct position in the array.
     * if(po.length>=2) {       //checks functions with at least x^1 (the derivative of a constant term is 0)
     * ans = new double[po.lengh-1]         //change the length of ans to po length -1
     * for(int i = po.lengh-1;i>=0;i--){    //go through all the values in po from the end
     * ans[i-1]=po[i]*i }                    //chanfe the values in ans to the power * the coefficient
     * }                                    //ends for block
     * return ans                           //return final derivative array
	 */
	public static double[] derivative (double[] po) {
		double [] ans = ZERO;
        if (po.length>=2) {
            ans = new double[po.length-1];
            for(int i=po.length-1;i>0;i--) {
                ans[i-1]=po[i]*i;
            }
        }
		return ans;
	}

    /**
     * this function gets 2 double arrays and return the longer array
     * @param p1 a double array
     * @param p2 a double array
     * @return the longer array, if they have the same length, return the first one (p1)
     * this function checks if p1.length is bigger or equals to p2.length and return p1 if so, else return p2.
     * double max_p [];
     * if(p1.length>=p2.length) {max_p=p1}
     * else {max_p = 2}
     * return max_p
     */
    private static double [] maxArray(double[] p1, double[] p2) {
        double [] max_p;
        if (p1.length>=p2.length) {
            max_p = p1;
        }
        else  {
            max_p = p2; }
        return max_p;
        }
    /**
     * this function gets 2 double arrays and return the shorter array
     * @param p1 a double array
     * @param p2 a double array
     * @return the shorter array, if they have the same length, return the second one (p2)
     * this function checks if p1.length is less than p2.length and return p1 if so, else return p2.
     * double min_p [];
     * if(p1.length<p2.length) {min_p=p1}
     * else {min_p = 2}
     * return min_p
     */
    private static double [] minArray(double[] p1, double[] p2) {
        double [] min_p;
        if (p1.length<p2.length) {
            min_p = p1;
        }
        else  {
            min_p = p2; }
        return min_p;
    }

    /**
     * this function checks if all the array's values are zero
     * @param p double array
     * @return false if at least one of the array's values is not 0, else return true.
     * this function check for each value of the array if it's not 0,if finds one, returns true
     * for (int i = 0; i < p.length; i++) {
     *      if (p[i] != 0) {
     *          return false;}
     * }
     * return true;
     */
    private static boolean isAllZero(double[] p) {
        for (int i = 0; i < p.length; i++) {
            if (p[i] != 0) {
                return false;
            }
        }
        return true;
    }

    /**
     * this function compute the equation of a straight line using the Lagrange Formula (detailed in the link above)
     * @param xx an array with values of x1,x2
     * @param yy an array with values of y1,y2
     * @return an array representing the equation of the line that connects the 3 points
     * this function computes the coefficients of the polynom and inserts them to the array
     * double denom = (xx[0] - xx[1]) * (xx[0] - xx[2]) * (xx[1] - xx[2]); //computing the denominator for finding the equation
     * if (denom == 0) {return null;} //making sure that the denominator is not 0
     *  double A = (xx[2] * (yy[1] - yy[0]) + xx[1] * (yy[0] - yy[2]) + xx[0] * (yy[2] - yy[1])) / denom;
     *  double B = (xx[2]*xx[2] * (yy[0] - yy[1]) + xx[1]*xx[1] * (yy[2] - yy[0]) + xx[0]*xx[0] * (yy[1] - yy[2])) / denom;
     *  double C = (xx[1] * xx[2] * (xx[1] - xx[2]) * yy[0] + xx[2] * xx[0] * (xx[2] - xx[0]) * yy[1] + xx[0] * xx[1] * (xx[0] - xx[1]) * yy[2]) / denom;
     *  double [] ans = {C,B,A}; //inserting the final coefficients into the array
     *  return ans
     */
    private static double [] polyForThree(double[] xx, double[] yy) {
        double denom = (xx[0] - xx[1]) * (xx[0] - xx[2]) * (xx[1] - xx[2]);
        if (denom == 0) {return null;}
        double A = (xx[2] * (yy[1] - yy[0]) + xx[1] * (yy[0] - yy[2]) + xx[0] * (yy[2] - yy[1])) / denom;
        double B = (xx[2]*xx[2] * (yy[0] - yy[1]) + xx[1]*xx[1] * (yy[2] - yy[0]) + xx[0]*xx[0] * (yy[1] - yy[2])) / denom;
        double C = (xx[1] * xx[2] * (xx[1] - xx[2]) * yy[0] + xx[2] * xx[0] * (xx[2] - xx[0]) * yy[1] + xx[0] * xx[1] * (xx[0] - xx[1]) * yy[2]) / denom;
        double [] ans = {C,B,A};
        return ans;
    }

    /**
     * this function compute the equation of a straight line using the Slope-Intercept Form
     * @param xx an array with values of x1,x2
     * @param yy an array with values of y1,y2
     * @return an array representing the equation of the line that connects both points
     * this function computes the slope - m and the intercept - b and inserts them to the final array.
     * double numer = (yy[1]-yy[0])        //computing the numerator for finding the gradient
     * double denom = (xx[1] - xx[0])     //computing the denominator for finding the Gradient
     * if (denom == 0) {return null}       //making sure that the denominator is not 0
     * double m = numer/denom              //computing the Gradient
     * double b = yy[1]-m*xx[1]            //computing the intercept
     * double [] ans = {b,m}               //insert to final array, when m is the
     * return ans
     */
    private static double [] polyForTwo(double[] xx, double[] yy) {
        double numer = (yy[1]-yy[0]);
        double denom = (xx[1] - xx[0]);
        if (denom == 0) {return null;}
        double m = numer/denom;
        double b = yy[1]-m*xx[1];
        double [] ans = {b,m};
        return ans;
    }

    /**
     * this function creates a double array with fixed values 0-size ({0.0,1.0,2.0,...,size.0})
     * @param size an integer setting the size of the array
     * @return a double array with fixed values
     * this function creates a double array and inserts i in arr[i] (which converts it to double)
     * double [] arr = new double [size];
     * for (int i = 0; i < size; i++) {arr[i] = i;}
     * return arr
     */
    private static double [] createFinalArr (int size) {
        double [] arr = new double [size];
        for (int i = 0; i < size; i++) {
            arr[i] = i;
        }
        return arr;
    }
}

