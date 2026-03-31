import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { GLTFLoader } from 'three/addons/loaders/GLTFLoader.js';
import Stats from 'three/addons/libs/stats.module.js';
import GUI from 'lil-gui';
import createModule from './fftocean_wasm.js';

let scene, camera, renderer;
let geometry, cube, material;
let coordinateArray;
let xyz_geometry;
let mesh;       // instanced mesh
let stats;
let ocean;
const L = 300;
const N = 256;

const ocean_settings = {
    height_scale: 2000.0,
    choppy_coefficient: 0.0,
}


createModule().then((Module) => {

    const setup_params = {
        wx: 10.0,
        wz: 10.0,
        initialize: function() {
            if (ocean instanceof Module.FFTOcean) {ocean.delete();}
            ocean = new Module.FFTOcean(
                L, N, setup_params.wx, setup_params.wz, 
                0.005, ocean_settings.choppy_coefficient, ocean_settings.height_scale);
        }
    };
    
    ocean = new Module.FFTOcean(L, N, 10.0, 10.0, 0.005, 1.0, 2000);
    setup_three();
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
        const axesHelper = new THREE.AxesHelper(10);
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
        folder.add(ocean_settings, 'choppy_coefficient', 0, 1000, 100).onChange(value => {ocean.set_choppy_coefficient(value)});

        const setup_param_folder = gui.addFolder('Setup Params');
        setup_param_folder.add(setup_params, 'wx', 0, 30);
        setup_param_folder.add(setup_params, 'wz', 0, 30);
        setup_param_folder.add(setup_params, 'initialize');
    }
    function setup_three() {
        scene = new THREE.Scene();
        camera = new THREE.PerspectiveCamera(
            75, window.innerWidth / window.innerHeight, 0.1, 1000
        );
        camera.position.set(5, 5, 5);
        camera.lookAt(0, 0, 0);

        renderer = new THREE.WebGLRenderer({antialias: true});
        renderer.setSize(window.innerWidth, window.innerHeight);
        renderer.toneMapping = THREE.ACESFilmicToneMapping;
        renderer.toneMappingExposure = 1.0; // 露出調整。暗ければ 1.5 くらいに。
        document.body.appendChild(renderer.domElement);

        geometry = new THREE.PlaneGeometry(L, L, N-1, N-1);
        material = new THREE.MeshNormalMaterial({wireframe: true, side: THREE.DoubleSide});

        //material = new THREE.MeshStandardMaterial({
        //    color: 0x001eff,       // ベースとなる水の色（深い青緑）
        //    metalness: 1.0,        // 金属光沢（鏡のような反射にするため高めに）
        //    roughness: 0.1,        // 表面の粗さ（滑らかにするため低めに）
        //    wireframe: false,       // ワイヤーフレームをオフにする（形を見るならtrueでもOK）
        //    side: THREE.DoubleSide // 両面を描画
        //});
        const oceanMesh = new THREE.Mesh(geometry, material);
        scene.add(oceanMesh);

        // 太陽光（平行光源）を追加
        const sunLight = new THREE.DirectionalLight(0xffffff, 1.0); // 白色の光、強さ1.0
        sunLight.position.set(100, 100, 100); // 光の来る方向
        scene.add(sunLight);

        // 環境光（影の部分が真っ黒になるのを防ぐ）を追加
        const ambientLight = new THREE.AmbientLight(0xffffff, 2.0); 
        scene.add(ambientLight);
        const controls = new OrbitControls(camera, renderer.domElement);
    }

    function animate() {
        requestAnimationFrame(animate);

        const time = performance.now() / 1000;
        ocean.Update(time);

        const height_ptr = ocean.get_xyz_ptr();
        const grad_ptr = ocean.get_gxyz_ptr();
        const heightArray = new Float32Array(Module.HEAPF32.buffer, height_ptr, N*N*3);
        const gradArray = new Float32Array(Module.HEAPF32.buffer, grad_ptr, N*N*3);

        //const posAttr = geometry.attributes.position;
        geometry.attributes.position.array.set(heightArray);
        geometry.attributes.normal.array.set(gradArray);
        geometry.attributes.position.needsUpdate = true;
        geometry.attributes.normal.needsUpdate = true;
        stats.update();
        renderer.render(scene, camera);
    }


});