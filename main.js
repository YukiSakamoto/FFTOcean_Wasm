import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { GLTFLoader } from 'three/addons/loaders/GLTFLoader.js';
import Stats from 'three/addons/libs/stats.module.js';
import GUI from 'lil-gui';
import createModule from './fftocean_wasm.js';

let scene, camera, renderer;
let geometry, cube, material;
let axesHelper;
let prev_time_virtual = 0;
let prev_time_real = 0;
let coordinateArray;
let xyz_geometry;
let mesh;       // instanced mesh
let stats;
let ocean;
const L = 300;
const N = 256;
//const L = 1000;
//const N = 512;

const ocean_settings = {
    height_scale: 2000.0,
    choppy_coefficient: 0.0,

    speed_factor: 1.0,
    show_axes: true,
}


createModule().then((Module) => {

    const setup_params = {
        wx: 10.0,
        wz: 10.0,
        A: 0.005,
        initialize: init_FFTOcean
    };

    function init_FFTOcean() {
        if (ocean instanceof Module.FFTOcean) {ocean.delete();}
        ocean = new Module.FFTOcean(
            L, N, setup_params.wx, setup_params.wz, 
            setup_params.A, ocean_settings.choppy_coefficient, ocean_settings.height_scale);
        prev_time_real = 0.0;
        prev_time_virtual = 0.0;

        const height_ptr = ocean.get_xyz_ptr();
        const grad_ptr = ocean.get_gxyz_ptr();
        const heightArray = new Float32Array(Module.HEAPF32.buffer, height_ptr, N*N*3);
        const gradArray = new Float32Array(Module.HEAPF32.buffer, grad_ptr, N*N*3);
        geometry.setAttribute('position', new THREE.BufferAttribute(heightArray, 3));
        geometry.setAttribute('normal', new THREE.BufferAttribute(gradArray, 3));
        geometry.attributes.position.setUsage(THREE.DynamicDrawUsage);
        geometry.attributes.normal.setUsage(THREE.DynamicDrawUsage);
        geometry.attributes.position.needsUpdate = true;
        geometry.attributes.normal.needsUpdate = true;
    }
    
    //ocean = new Module.FFTOcean(L, N, 10.0, 10.0, 0.005, 1.0, 2000);
    setup_three();
    init_FFTOcean();
    setup_stats();
    setup_axisHelper();
    setup_lilgui();
    animate();
    window.addEventListener('resize', () => {
        camera.aspect = window.innerWidth / window.innerHeight;
        camera.updateProjectionMatrix();
        renderer.setSize(window.innerWidth, window.innerHeight);
    })
    function setup_axisHelper() {
        axesHelper = new THREE.AxesHelper(10);
        scene.add(axesHelper);
    }
    function setup_stats() {
        stats = new Stats();
        stats.showPanel(0);
        document.body.appendChild(stats.dom);
    }
    function setup_lilgui() {
        const gui = new GUI();
        const folder = gui.addFolder('settings');
        folder.add(ocean_settings, 'height_scale', 0, 5000, 1).onChange(value => {ocean.set_height_scale(value)});
        folder.add(ocean_settings, 'choppy_coefficient', 0, 1000, 1).onChange(value => {ocean.set_choppy_coefficient(value)});
        folder.add(ocean_settings, 'speed_factor', 0, 10, 1);
        folder.add(ocean_settings, 'show_axes').onChange(value => {axesHelper.visible = value})

        const setup_param_folder = gui.addFolder('Setup Params');
        setup_param_folder.add(setup_params, 'wx', 0, 30);
        setup_param_folder.add(setup_params, 'wz', 0, 30);
        setup_param_folder.add(setup_params, 'A', 0, 0.1, 0.005);
        setup_param_folder.add(setup_params, 'initialize');
    }
    function setup_three() {
        scene = new THREE.Scene();
        camera = new THREE.PerspectiveCamera(
            75, window.innerWidth / window.innerHeight, 0.1, 1000
        );
        camera.position.set(0, 50, -100);
        camera.lookAt(0, 0, 0);

        renderer = new THREE.WebGLRenderer({antialias: true});
        renderer.setSize(window.innerWidth, window.innerHeight);
        renderer.toneMapping = THREE.ACESFilmicToneMapping;
        renderer.toneMappingExposure = 1.0; // 露出調整。暗ければ 1.5 くらいに。
        document.body.appendChild(renderer.domElement);

        //geometry = new THREE.PlaneGeometry(L, L, N-1, N-1);
        //material = new THREE.MeshNormalMaterial({wireframe: true, side: THREE.DoubleSide});

        const indices = [];
        for (let i = 0; i < N - 1; i++) {
            for (let j = 0; j < N - 1; j++) {
                const a = i * N + j;
                const b = (i + 1) * N + j;
                const c = i * N + (j + 1);
                const d = (i + 1) * N + (j + 1);

                indices.push(a, c, b);
                indices.push(b, c, d);
            }
        }
        geometry = new THREE.BufferGeometry();
        geometry.setIndex(indices);
        //geometry.setAttribute('position', new THREE.BufferAttribute(heightArray, 3));
        //geometry.setAttribute('normal', new THREE.BufferAttribute(gradArray, 3));

        material = new THREE.MeshStandardMaterial({
            color: 0x00bfff,       // ベースとなる水の色（深い青緑）
            metalness: 0.9,        // 金属光沢（鏡のような反射にするため高めに）
            roughness: 0.5,        // 表面の粗さ（滑らかにするため低めに）
            //wireframe: true,       // ワイヤーフレームをオフにする（形を見るならtrueでもOK）
            side: THREE.FrontSide
        });
        const oceanMesh = new THREE.Mesh(geometry, material);
        scene.add(oceanMesh);

        // 太陽光（平行光源）を追加
        const sunLight = new THREE.DirectionalLight(0xffffff, 1.0); // 白色の光、強さ1.0
        sunLight.position.set(0, 100, 100); // 光の来る方向
        scene.add(sunLight);

        // 環境光（影の部分が真っ黒になるのを防ぐ）を追加
        const ambientLight = new THREE.AmbientLight(0xffffff, 2.0); 
        scene.add(ambientLight);
        const controls = new OrbitControls(camera, renderer.domElement);
        scene.environment = scene.background;
    }

    function animate() {
        requestAnimationFrame(animate);

        const time_real = performance.now();
        const dt_real = time_real - prev_time_real;
        const time_virtual = prev_time_virtual + dt_real / 1000 * ocean_settings.speed_factor;
        ocean.Update(time_virtual);
        // store time as previous values;
        prev_time_virtual = time_virtual;
        prev_time_real = time_real;

        geometry.attributes.position.needsUpdate = true;
        geometry.attributes.normal.needsUpdate = true;
        stats.update();
        renderer.render(scene, camera);
    }


});