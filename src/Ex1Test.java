import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 *  * Introduction to Computer Science 2026, Ariel University,
 *  * Ex1: arrays, static functions and JUnit
 *
 * This JUnit class represents a JUnit (unit testing) for Ex1-
 * It contains few testing functions for the polynomial functions as define in Ex1.
 * Note: you should add additional JUnit testing functions to this class.
 *
 * @author boaz.ben-moshe
 */

class Ex1Test {
	static final double[] P1 ={2,0,3, -1,0}, P2 = {0.1,0,1, 0.1,3};
	static double[] po1 = {2,2}, po2 = {-3, 0.61, 0.2};;
	static double[] po3 = {2,1,-0.7, -0.02,0.02};
	static double[] po4 = {-3, 0.61, 0.2};
	
 	@Test
	/**
	 * Tests that f(x) == poly(x).
	 */
	void testF() {
		double fx0 = Ex1.f(po1, 0);
		double fx1 = Ex1.f(po1, 1);
		double fx2 = Ex1.f(po1, 2);
		assertEquals(fx0, 2, Ex1.EPS);
		assertEquals(fx1, 4, Ex1.EPS);
		assertEquals(fx2, 6, Ex1.EPS);
	}
	@Test
	/**
	 * Tests that p1(x) + p2(x) == (p1+p2)(x)
	 */
	void testF2() {
		double x = Math.PI;
		double[] po12 = Ex1.add(po1, po2);
		double f1x = Ex1.f(po1, x);
		double f2x = Ex1.f(po2, x);
		double f12x = Ex1.f(po12, x);
		assertEquals(f1x + f2x, f12x, Ex1.EPS);
	}
    @Test
    void PolynomFromTreePoints(){
        double[] xx1 = {1,3,5};
        double[] yy1 = {2,4,6};
        double [] polyOfXxYy1 = Ex1.PolynomFromPoints(xx1,yy1);
        double [] expected_1 = {1,1,0};
        assertArrayEquals(polyOfXxYy1,expected_1,Ex1.EPS);
        double[] xx2 = {14,3,0};
        double[] yy2 = {61,6,33};
        double [] polyOfXxYy2 = Ex1.PolynomFromPoints(xx2,yy2);
        double [] expected_2 = {33,-12,1};
        assertArrayEquals(polyOfXxYy2,expected_2,Ex1.EPS);
    }
    @Test
    void PolynomFromTwoPoints() {
        double[] xx1 = {1,0};
        double[] yy1 = {2,6};
        double [] polyOfXxYy1 = Ex1.PolynomFromPoints(xx1,yy1);
        double [] expected_1 = {6,-4};
        assertArrayEquals(polyOfXxYy1,expected_1,Ex1.EPS);
        double [] xx2={7,8};
        double[] yy2 = {2,0};
        double [] polyOfXxYy2 = Ex1.PolynomFromPoints(xx2, yy2);
        double [] expected_2 = {16,-2};
        assertArrayEquals(polyOfXxYy2, expected_2,Ex1.EPS);
    }
    @Test
    void testAddSimple(){
        double [] p1 = {1,2,3};
        double [] p2 = {1,2,-3};
        double [] sum_p12 = Ex1.add(p1, p2);
        double [] expected_sum12 = {2,4,0};
        assertArrayEquals(expected_sum12,sum_p12,Ex1.EPS);
        double [] p3 = {5,6.5,23};             //23x^2+6.5x+5
        double [] p4 = {3.4,-2,5,4};           //4x^3+5x^2-2x+3.4
        double [] sum_p34 = Ex1.add(p3,p4);
        double [] expected_sum34 = {8.4,4.5,28,4};        //4x^3+28x^2+4.5x+8.4
        assertArrayEquals(sum_p34, expected_sum34,Ex1.EPS);
        double [] p5 = {0};                 //0
        double [] big_p6 = {1000,-10000,100000};                 //100000x^2-10000x+1000
        double [] sum_p56 = Ex1.add(p5,big_p6);
        double [] expected_sum_p56 = {1000,-10000,100000};       //p4+0=p4
        assertArrayEquals(sum_p56, expected_sum_p56,Ex1.EPS);
    }
	@Test
	/**
	 * Tests that p1+p2+ (-1*p2) == p1
	 */
	void testAdd() {
		double[] p12 = Ex1.add(po1, po2);
		double[] minus1 = {-1};
		double[] pp2 = Ex1.mul(po2, minus1);
		double[] p1 = Ex1.add(p12, pp2);
		assertTrue(Ex1.equals(p1, po1));
	}
	@Test
	/**
	 * Tests that p1+p2 == p2+p1
	 */
	void testAdd2() {
		double[] p12 = Ex1.add(po1, po2);
		double[] p21 = Ex1.add(po2, po1);
		assertTrue(Ex1.equals(p12, p21));
	}
	@Test
	/**
	 * Tests that p1+0 == p1
	 */
	void testAdd3() {
		double[] p1 = Ex1.add(po1, Ex1.ZERO);
		assertTrue(Ex1.equals(p1, po1));
	}
    @Test
    void testMulSimple() {
        double[] p1 = {1, 2, 3};
        double[] p2 = {1, 2};
        double[] mul_p12 = Ex1.mul(p1, p2);
        double[] expected_mul12 = {1, 4, 7, 6};
        assertArrayEquals(expected_mul12, mul_p12, Ex1.EPS);
    }
    @Test
    void testComplex() {
        double [] p3 = {10,2,2.2,-1};             //-x^3+2.2x^2+2x+10
        double [] p4 = {3.9,-2.5,0,11};           //11x^3-2.5x+3.9
        double [] mul_p34 = Ex1.mul(p3,p4);
        double [] expected_mul34 = {39,-17.2,3.58,100.6,24.5,24.2,-11};        //-11x^6+24.2x^5+24.5x^4+100.6x^3+3.58x^2-17.2x+39
        assertArrayEquals(mul_p34, expected_mul34,Ex1.EPS);

    }@Test
    void testbig2() {
        double [] p1 = {379,4493,539}; //539x^2+4493x+379
        double [] p2 = {111,222,333};  //333x^2+222x+111
        double [] big_p12 = Ex1.mul(p1,p2);
        double [] expected_mul12 = {42069,582861,1183482,1615827,179487};
        assertArrayEquals(big_p12, expected_mul12,Ex1.EPS);
        double [] p5 = {0};                 //0
        double [] big_p6 = {1000,-10000,100000};                 //100000x^2-10000x+1000
        double [] mul_p56 = Ex1.mul(p5,big_p6);
        double [] expected_mul_p56 = Ex1.ZERO;       //p4*0=p4
        assertArrayEquals(mul_p56, expected_mul_p56,Ex1.EPS);
    }
	@Test
	/**
	 * Tests that p1*0 == 0
	 */
	void testMul1() {
		double[] p1 = Ex1.mul(po1, Ex1.ZERO);
		assertTrue(Ex1.equals(p1, Ex1.ZERO));
	}
	@Test
	/**
	 * Tests that p1*p2 == p2*p1
	 */
	void testMul2() {
		double[] p12 = Ex1.mul(po1, po2);
		double[] p21 = Ex1.mul(po2, po1);
		assertTrue(Ex1.equals(p12, p21));
	}
	@Test
	/**
	 * Tests that p1(x) * p2(x) = (p1*p2)(x),
	 */
	void testMulDoubleArrayDoubleArray() {
		double[] xx = {0,1,2,3,4.1,-15.2222};
		double[] p12 = Ex1.mul(po1, po2);
		for(int i = 0;i<xx.length;i=i+1) {
			double x = xx[i];
			double f1x = Ex1.f(po1, x);
			double f2x = Ex1.f(po2, x);
			double f12x = Ex1.f(p12, x);
			assertEquals(f12x, f1x*f2x, Ex1.EPS);
		}
	}
    @Test
    void testDerivativeSimple() {
        double [] p1 = {1,2,3};             //3x^2x+1
        double [] dp1 = Ex1.derivative(p1);
        double [] expected_dp1 = {2,6};        //6x+2
        assertArrayEquals(dp1,expected_dp1,Ex1.EPS);
        double [] p2 = {3.4,-2.7,5,4};      //4x^3+7.5x^2-2.7x+3.4
        double [] dp2 = Ex1.derivative(p2);
        double [] expected_dp2 = {-2.7,10,12};  //12x^2+10x-2.7
        assertArrayEquals(dp2,expected_dp2,Ex1.EPS);
        double [] p3 = {0};                 //0
        double [] dp3 = Ex1.derivative(p3);
        double [] expected_dp3 = {0};       //0
        assertArrayEquals(dp3,expected_dp3,Ex1.EPS);
        double [] big_p4 = {1000,-10000,100000};                 //100000x^2-10000x+1000
        double [] big_dp4 = Ex1.derivative(big_p4);
        double [] expected_big_dp4 = {-10000,200000};       //200000x-10000
        assertArrayEquals(big_dp4, expected_big_dp4,Ex1.EPS);
    }
	@Test
	/**
	 * Tests a simple derivative examples - till ZERO.
	 */
	void testDerivativeArrayDoubleArray() {
		double[] p = {1,2,3}; // 3X^2+2x+1
		double[] pt = {2,6}; // 6x+2
		double[] dp1 = Ex1.derivative(p); // 2x + 6
		double[] dp2 = Ex1.derivative(dp1); // 2
		double[] dp3 = Ex1.derivative(dp2); // 0
		double[] dp4 = Ex1.derivative(dp3); // 0
		assertTrue(Ex1.equals(dp1, pt));
		assertTrue(Ex1.equals(Ex1.ZERO, dp3));
		assertTrue(Ex1.equals(dp4, dp3));
	}
	@Test
	/** 
	 * Tests the parsing of a polynom in a String like form.
	 */
	public void testFromString() {
		double[] p = {-1.1,2.3,-3.1,-7}; // 3.1X^2+ 2.3x -1.1
		String sp2 = "-7.0x^3-3.1x^2+2.3x-1.1";
		String sp = Ex1.poly(p);
		double[] p1 = Ex1.getPolynomFromString(sp);
		double[] p2 = Ex1.getPolynomFromString(sp2);
		boolean isSame1 = Ex1.equals(p1, p);
		boolean isSame2 = Ex1.equals(p2, p);
		if(!isSame1) {fail();}
		if(!isSame2) {fail();}
		assertEquals(sp, Ex1.poly(p1));
	}
    @Test
    void testFromStringSimple() {
        double[] p = {-1.1,2.3,3.1}; // 3.1X^2+ 2.3x -1.1
        String e_sp = "3.1x^2+2.3x-1.1";
        assertEquals(e_sp,Ex1.poly(p));
        double[] p2 = {-5,0,5.7,3.1,0,-5,3}; // 3x^6+-5x^5+3.1x^3+5.7x^2-5
        String e_sp2 = "3.0x^6-5.0x^5+3.1x^3+5.7x^2-5.0";
        String sp2 = Ex1.poly(p2);
        assertEquals(e_sp2, sp2);
    }

	@Test
	/**
	 * Tests the equality of pairs of arrays.
	 */
	public void testEquals() {
		double[][] d1 = {{0}, {1}, {1,2,0,0}};
		double[][] d2 = {Ex1.ZERO, {1+ Ex1.EPS/2}, {1,2}};
		double[][] xx = {{-2* Ex1.EPS}, {1+ Ex1.EPS*1.2}, {1,2, Ex1.EPS/2}};
		for(int i=0;i<d1.length;i=i+1) {
			assertTrue(Ex1.equals(d1[i], d2[i]));
		}
		for(int i=0;i<d1.length;i=i+1) {
			assertFalse(Ex1.equals(d1[i], xx[i]));
		}
	}
    @Test
    void testLength () {
        double[]p = {0,1};
        double x1 = 0;
        double x2 = 3;
        int NoS = 3;
        double expected = 3 * Math.sqrt(2);
        assertEquals(expected, Ex1.length(p, x1, x2, NoS),Ex1.EPS);
        double[] p2 = {1, 3, -1, 2};
        double x3 = 0;
        double x4 = 2;
        int NoS2 = 1000;
        double expected2 = 18.1729;
        assertEquals(expected2, Ex1.length(p2, x3, x4, NoS2), Ex1.EPS);
    }


    @Test
    void testEqualsSimple() {
        double[] p1 = {5.0, 1.0, 0.0};
        double[] p2 = {5.0, 1.0};
        assertTrue(Ex1.equals(p1, p2));
        double[] p3 = {10.0};
        double[] p4 = {10.0 + Ex1.EPS / 3.0};
        assertTrue(Ex1.equals(p3, p4));
        double[] p5 = {-0.5, 4.25, 1.0, 0.0};
        double[] p6 = {-0.5, 4.25, 1.0};
        assertTrue(Ex1.equals(p5, p6));
    }
    @Test
    void testSameValueSimple() {
        double[] p1 = {1.0, 1.0};
        double[] p2 = {2.0, -1.0};
        double x1 = 0.0;
        double x2 = 5.0;
        double expectedX = 0.5;
        assertTrue(Ex1.sameValue(p1,p2,x1,x2,Ex1.EPS)<=0.51 && Ex1.sameValue(p1,p2,x1,x2,Ex1.EPS)>=0.4999);
    }

	@Test
	/**
	 * Tests is the sameValue function is symmetric.
	 */
	public void testSameValue2() {
		double x1=-4, x2=0;
		double rs1 = Ex1.sameValue(po1,po2, x1, x2, Ex1.EPS);
		double rs2 = Ex1.sameValue(po2,po1, x1, x2, Ex1.EPS);
		assertEquals(rs1,rs2, Ex1.EPS);
	}
	@Test
	/**
	 * Test the area function - it should be symmetric.
	 */
	public void testArea() {
		double x1=-4, x2=0;
		double a1 = Ex1.area(po1, po2, x1, x2, 100);
		double a2 = Ex1.area(po2, po1, x1, x2, 100);
		assertEquals(a1,a2, Ex1.EPS);
}
	@Test
	/**
	 * Test the area f1(x)=0, f2(x)=x;
     * This test fails because of the approximation error comes from the trapezoidal method (2.5<2.501953125)
	 */
	public void testArea2() {
		double[] po_a = Ex1.ZERO;
		double[] po_b = {0,1};
		double x1 = -1;
		double x2 = 2;
		double a1 = Ex1.area(po_a,po_b, x1, x2, 1);
		double a2 = Ex1.area(po_a,po_b, x1, x2, 2);
		double a3 = Ex1.area(po_a,po_b, x1, x2, 3);
		double a100 = Ex1.area(po_a,po_b, x1, x2, 100);
		double area =2.5;
		assertEquals(a1,area, Ex1.EPS);
		assertEquals(a2,area, Ex1.EPS);
		assertEquals(a3,area, Ex1.EPS);
		assertEquals(a100,area, Ex1.EPS);
	}
	@Test
	/**
	 * Test the area function.
     * This test fails because of the approximation error comes from the trapezoidal method (2.5<2.501953125)
	 */
	public void testArea3() {
		double[] po_a = {2,1,-0.7, -0.02,0.02};
		double[] po_b = {6, 0.1, -0.2};
		double x1 = Ex1.sameValue(po_a,po_b, -10,-5, Ex1.EPS);
		double a1 = Ex1.area(po_a,po_b, x1, 6, 8);
		double area = 58.5658;
		assertEquals(a1,area, Ex1.EPS);
	}
}
